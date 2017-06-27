

#include "headersBasic.h"

#include "StandardFEM.h"
#include "MathBasic.h"
#include "DataBlockTemplate.h"
#include "ComputerTime.h"
#include "MpapTime.h"
#include "TimeFunction.h"
#include "PlotVTK.h"
#include "PropertyTypeEnum.h"

#include "Global.h"
#include "MyString.h"
#include "FunctionsEssGrp.h"

#include "UnixGlobal.h"
#include <omp.h>
#include "petsc.h"


extern DomainTree         domain;
extern List<TimeFunction> timeFunction;
extern MpapTime           mpapTime;
extern ComputerTime       computerTime;
extern PlotVTK  plotvtk;


using namespace std;



//
int main(int argc, char* argv[])
{
  PetscErrorCode ierr;

  //PetscInitialize(NULL,NULL,(char *)0,NULL);

  PetscInitialize(NULL, NULL, "petsc_options.dat", NULL);

  // FEM using petsc library
  // to use parallel solvers

    int parm[8], ii, jj, iter, solv=8;
    parm[0] = 1;

    int  resln[]={1,1,1};
    int tis = 0 ;
    double  tol = 1.0e-5;
    double  rho = 0.8;

    vector<double>  stagParams(3);

    stagParams[0] = 0;
    stagParams[1] = 2;
    //stagParams[2] = 1.2/(1.0+1.2)/200.0;
    stagParams[2] = 0.1;

    files.projDir = "../src/parallel/" ;
    files.Ifile = "Ioutput";
    files.Ofile.free().append(files.Ifile)[0] = 'O';
    files.Tfile.free().append(files.Ifile)[0] = 'T';
    files.Pfile.free().append(files.Ifile)[0] = 'P';

    cout << "       project directory   : " << files.projDir << "\n";
    cout << "       input file name     : " << files.Ifile << "\n";
    cout << "       output file name    : " << files.Ofile << "\n";
    cout << "       time plot file name : " << files.Tfile << "\n";
    cout << "       eps plot file name  : " << files.Pfile << "\n\n";  


    // ==================================================
    //
    // prepare input data for the solid problem
    //
    // ==================================================

    StandardFEM  femObject;

    //std::ifstream  infile("../src/parallel/ILaplace2DEx1tria20");
    std::ifstream  infile("../src/parallel/ICylinderTriasMesh2");

    if(infile.fail())
    {
       cout << " Could not open the input file" << endl;
      exit(1);
    }
    //cout << " llllllllllll " << endl;
    femObject.SetDimension(2);

    femObject.setNdof(3);

    femObject.readFile( infile );
    //cout << " llllllllllll " << endl;

    femObject.SolnData.SetStaggeredParams(stagParams);

    femObject.SetVTKfilename("ILaplace2DEx1tria20");

    femObject.prepareInputData();
  
    cout << " Solid model preparation successful " << endl;

    femObject.setSolver(solv, parm);

    cout << " Solid solver preparation successful " << endl;
  


    // ==================================================
    //
    // set the time function(s)
    //
    // ==================================================


    timeFunction.add(new TimeFunction);

    timeFunction[0].fct.add(new TimeFunctionCore);

    timeFunction[0].fct[0].t0 = 0.0;
    timeFunction[0].fct[0].t1 = 10.0;

    timeFunction[0].fct[0].p[0] = 1.0;
    timeFunction[0].fct[0].p[1] = 0.0;
    timeFunction[0].fct[0].p[2] = 0.0;
    timeFunction[0].fct[0].p[3] = 0.0;
    timeFunction[0].fct[0].p[4] = 0.0;
    timeFunction[0].fct[0].p[5] = 0.0;

    timeFunction[0].fct[0].tp = 1.0e+30;

    for(ii=0; ii<timeFunction.n; ii++)
      timeFunction[ii].update();


    double  dt=1.0, tf=1.0, tCur;



    mpapTime.dtOK     = true;
    mpapTime.dt       = dt;
    mpapTime.dtMax    = dt;
    mpapTime.stack.free();
    mpapTime.stack.append(mpapTime.dt);
    if(mpapTime.dtMax < 1.e-15)
      mpapTime.dtMax = mpapTime.dt;

    tCur  = dt;

    femObject.SolnData.SetStaggeredParams(stagParams);

    int niter=2;

    femObject.setTimeParam();

    femObject.timeUpdate();

    for(iter=0; iter<niter; iter++)
    {
      //cout << " aaaaaaaaaaa " << endl;
      femObject.calcStiffnessAndResidual();

      //cout << " bbbbbbbbbbb " << endl;
      if( femObject.converged() )
        break;

      //cout << " ccccccccccc " << endl;
      femObject.factoriseSolveAndUpdate();

      //cout << " ddddddddddd " << endl;
      femObject.updateIterStep();
    }

    femObject.postProcess(0, 0, 1, 0, 0.0, 1.0, resln);

    ierr = PetscFinalize();//CHKERRQ(ierr);

    cout << " \n\n\n Propgram successful ... \n\n\n " << endl;


  return 0;
}











