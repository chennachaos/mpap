

#include "headersBasic.h"

#include "StandardFEM.h"
#include "MathBasic.h"
#include "DataBlockTemplate.h"
#include "ComputerTime.h"
#include "MpapTime.h"
#include "TimeFunction.h"
#include "PropertyTypeEnum.h"
#include "Global.h"
#include "MyString.h"
#include "FunctionsEssGrp.h"
#include "UnixGlobal.h"
#include "petsc.h"
#include <omp.h>


extern DomainTree         domain;
extern List<TimeFunction> timeFunction;
extern MpapTime           mpapTime;
extern ComputerTime       computerTime;


using namespace std;



//
int main(int argc, char* argv[])
{
  PetscErrorCode ierr;

  //PetscInitialize(NULL,NULL,(char *)0,NULL);

  PetscInitialize(NULL, NULL, "petsc_options.dat", NULL);

  // FSI with standard elements
  // Example from the paper 
  // Stabilization of explicit coupling in fluid–structure interaction involving fluid incompressibility
  // by
  // Erik Burman a, Miguel A. Fernández
  //

    int parm[5], ii, jj, iter, niter=10, solv=1;
    parm[0] = 1;

    int  resln[]={1,1,1};
    int tis = 3 ;
    double  tol = 1.0e-5;
    double  rho = 0.8;

    vector<double>  stagParams(3);

    stagParams[0] = 0;
    stagParams[1] = 1;
    //stagParams[2] = 1.2/(1.0+1.2)/200.0;
    stagParams[2] = 0.1;

    files.projDir = "../src/test/" ;
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

    StandardFEM  sfemSolid;


    std::ifstream  infileSolid("../src/test/IChannelSolid");

    if(infileSolid.fail())
    {
       cout << " Could not open the 'SOLID' input file" << endl;
      exit(1);
    }

    sfemSolid.SetDimension(2);

    sfemSolid.SetNDOF(2);

    sfemSolid.readFile( infileSolid );

    sfemSolid.SetControl(tis, tol, rho);

    sfemSolid.SolnData.SetTimeIncrementType(tis);

    sfemSolid.SolnData.SetRho(rho);

    sfemSolid.SolnData.SetStaggeredParams(stagParams);
    
    sfemSolid.SetPhysicsTypetoSolid();

    sfemSolid.SetVTKfilename("FSIsolid");

    //sfemSolid.prepareInputData();
  
    cout << " Solid model preparation successful " << endl;

    sfemSolid.setSolver(solv, parm);

    cout << " Solid solver preparation successful " << endl;
  
    // ==================================================
    //
    // prepare input data for the fluid problem
    //
    // ==================================================

    StandardFEM  sfemFluid;

    sfemFluid.SetDimension(2);

    sfemFluid.SetNDOF(3);

    std::ifstream  infileFluid("../src/test/IChannelFluid");

    if(infileFluid.fail())
    {
       cout << " Could not open the 'FLUID' input file" << endl;
      exit(1);
    }

    sfemFluid.readFile( infileFluid );

    sfemFluid.SetControl(tis, tol, rho);

    sfemFluid.SolnData.SetTimeIncrementType(tis);

    sfemFluid.SolnData.SetRho(rho);

    sfemFluid.SetPhysicsTypetoFluid();

    sfemFluid.SetVTKfilename("FSIfluid");

    //sfemFluid.prepareInputData();
  
    cout << " Fluid model preparation successful " << endl;

    sfemFluid.setSolver(solv, parm);

    cout << " Fluid solver preparation successful " << endl;
  

    // ==================================================
    //
    // prepare input data for the mesh motion problem
    //
    // ==================================================

    StandardFEM  sfemMesh;

    sfemMesh.SetDimension(2);

    sfemMesh.SetNDOF(2);

    std::ifstream  infileMesh("../src/test/IChannelMesh");

    if(infileMesh.fail())
    {
      cout << " Could not open the 'FLUID' input file" << endl;
      exit(1);
    }

    sfemMesh.readFile( infileMesh );

    sfemMesh.SetControl(0, tol, rho);

    sfemMesh.SolnData.SetTimeIncrementType(0);

    sfemMesh.SolnData.SetRho(rho);

    sfemMesh.SetPhysicsTypetoSolid();

    sfemMesh.SetVTKfilename("FSIMesh");

    //sfemFluid.prepareInputData();
  
    cout << " Mesh model preparation successful " << endl;

    sfemMesh.setSolver(solv, parm);

    cout << " Mesh solver preparation successful " << endl;


    // ==================================================
    //
    // set the interface nodes
    //
    // ==================================================

    vector<int>  intfFluid(101), intfSolid(101);

    for(ii=0; ii<101; ii++)
    {
      intfFluid[ii] = 1010+ii;
      intfSolid[ii] = ii;
    }


    // ==================================================
    //
    // set the time function(s)
    //
    // ==================================================


    timeFunction.add(new TimeFunction);

    timeFunction[0].fct.add(new TimeFunctionCore);

    timeFunction[0].fct[0].t0 = 0.0;
    timeFunction[0].fct[0].t1 = 0.005;

    timeFunction[0].fct[0].p[0] = 1.0;
    timeFunction[0].fct[0].p[1] = 0.0;
    timeFunction[0].fct[0].p[2] = 0.0;
    timeFunction[0].fct[0].p[3] = 0.0;
    timeFunction[0].fct[0].p[4] = 0.0;
    timeFunction[0].fct[0].p[5] = 0.0;

    timeFunction[0].fct[0].tp = 1.0e+30;

    timeFunction[0].fct.add(new TimeFunctionCore);

    timeFunction[0].fct[1].t0 = 0.005;
    timeFunction[0].fct[1].t1 = 10.0;

    timeFunction[0].fct[1].p[0] = 0.0;
    timeFunction[0].fct[1].p[1] = 0.0;
    timeFunction[0].fct[1].p[2] = 0.0;
    timeFunction[0].fct[1].p[3] = 0.0;
    timeFunction[0].fct[1].p[4] = 0.0;
    timeFunction[0].fct[1].p[5] = 0.0;

    timeFunction[0].fct[1].tp = 1.0e+30;

    for(ii=0; ii<timeFunction.n; ii++)
      timeFunction[ii].update();


    double  dt=0.0001, tf=0.01, tCur;



    mpapTime.dtOK     = true;
    mpapTime.dt       = dt;
    mpapTime.dtMax    = dt;
    mpapTime.stack.free();
    mpapTime.stack.append(mpapTime.dt);
    if(mpapTime.dtMax < 1.e-15)
      mpapTime.dtMax = mpapTime.dt;

    niter = 10;
    tCur  = dt;

    sfemSolid.SolnData.SetStaggeredParams(stagParams);
    sfemFluid.SolnData.SetStaggeredParams(stagParams);

    while(tCur <= tf)
    {
        mpapTime.cur = tCur;

        // compute force predictor on the solid
        ///////////////////////////

        sfemSolid.setTimeParam();

        sfemSolid.timeUpdate();

        // solve solid problem
        ///////////////////////////

        cout << " solving solid problem ... " << endl;

        sfemSolid.SolveStep(niter);

        //printVector(sfemSolid.SolnData.var1Cur);

        cout << " solving solid problem ... DONE " << endl;

        sfemSolid.writeNodalData();

        sfemSolid.postProcess(0, 0, 1, 0, 0.0, 1.0, resln);


        //
        // mesh solver
        ///////////////////////////

        sfemMesh.SolnData.var1applied.setZero();
        for(ii=0; ii<intfFluid.size(); ii++)
        {
          for(jj=0; jj<2; jj++)
          {
            //cout << ii << '\t' << sfemSolid.SolnData.var1[intfSolid[ii]*2+jj] << endl;
            sfemMesh.SolnData.var1applied[intfFluid[ii]*2+jj] = sfemSolid.SolnData.var1[intfSolid[ii]*2+jj];
          }
        }

        sfemMesh.setTimeParam();

        sfemMesh.timeUpdate();

        cout << " solving mesh problem ... " << endl;

        sfemMesh.SolveStep(niter);

        cout << " solving mesh problem ... DONE " << endl;


        for(ii=0; ii<sfemMesh.GeomData.NodePosCur.size(); ii++)
        {
          sfemFluid.GeomData.NodePosCur[ii][0] = sfemMesh.GeomData.NodePosCur[ii][0];
          sfemFluid.GeomData.NodePosCur[ii][1] = sfemMesh.GeomData.NodePosCur[ii][1];

          sfemFluid.GeomData.NodePosNew[ii][0] = sfemMesh.GeomData.NodePosNew[ii][0];
          sfemFluid.GeomData.NodePosNew[ii][1] = sfemMesh.GeomData.NodePosNew[ii][1];
        }
        //
        // update fluid mesh and 
        // velocity BCs for the fluid problem
        ///////////////////////////

        sfemFluid.SolnData.var1applied.setZero();
        for(ii=0; ii<intfFluid.size(); ii++)
        {
          for(jj=0; jj<2; jj++)
          {
            //cout << ii << '\t' << sfemSolid.SolnData.var1DotCur[intfSolid[ii]*2+jj] << endl;
            sfemFluid.SolnData.var1[intfFluid[ii]*3+jj] = sfemSolid.SolnData.var1Dot[intfSolid[ii]*2+jj];
            //sfemFluid.SolnData.var1applied[intfFluid[ii]*3+jj] = sfemSolid.SolnData.var1Dot[intfSolid[ii]*2+jj];
          }
        }

        cout << " kkkkkkkkkkk " << endl;

        // solve fluid problem
        ///////////////////////////

        sfemFluid.setTimeParam();

        sfemFluid.timeUpdate();

        cout << " solving fluid problem ... " << endl;

        sfemFluid.SolveStep(niter);

        cout << " solving fluid problem ... DONE " << endl;

        //sfemFluid.writeNodalData();

        sfemFluid.postProcess(0, 0, 1, 0, 0.0, 1.0, resln);

        ///////////////////////////
        // compute the force corrector
        ///////////////////////////

        sfemSolid.SolnData.forceTemp.setZero();
        for(ii=0; ii<intfFluid.size(); ii++)
        {
          for(jj=0; jj<2; jj++)
          {
            sfemSolid.SolnData.forceTemp[intfSolid[ii]*2+jj] = sfemFluid.SolnData.reac[intfFluid[ii]*3+jj];
          }
        }

        sfemSolid.SolnData.interpolateForce();

        tCur += dt;
    }

    ierr = PetscFinalize();//CHKERRQ(ierr);

    cout << " \n\n\n Propgram successful ... \n\n\n " << endl;

    return 0;
}











