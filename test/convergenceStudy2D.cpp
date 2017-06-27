

#include "headersBasic.h"

#include "HBSplineFEM.h"
#include "HBSplineCutFEM.h"
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

#include "petsc.h"

#include <omp.h>
//#include "../../../../../../../usr/include/stdlib.h"


extern DomainTree         domain;
extern List<TimeFunction> timeFunction;
extern MpapTime           mpapTime;
extern ComputerTime       computerTime;
extern PlotVTK  plotvtk;


using namespace std;

/*
int main(int argc, char* argv[])
{
 // Advection-Diffusion 1D

  
  double  PEN;
  double  tol1 = 1.0e-6;
  double  rho1 = 0.8;
  double  NitscheFact;
  bool  isNitsche;
  
  int ne = 10, degree=1;

  PEN  = 1.0e3;
  PEN = 0.0;
  
  isNitsche = false;
  isNitsche = true;

  NitscheFact = -1.0;


  if(argc == 0)
    cerr << " Input data " << endl;
  else
  {
    ne  = atoi(argv[1]);
    degree  = atoi(argv[2]);
  }


  HBSplineFEM  hbs;

  //cout << "Hello, world! " << endl;
  
  hbs.SetDimension(1);

  hbs.SetOrigin(0.0, 0.0);
  
  hbs.SetGridDimensions(1.0, 1.0);
  
  hbs.SetBSplineDegree(degree, degree);
  
  hbs.SetNumberOfElements(ne, ne);


  hbs.AddDirichletBCs(1, 1, 0.0, PEN, isNitsche, NitscheFact);
  hbs.AddDirichletBCs(2, 1, 0.0, PEN, isNitsche, NitscheFact);
  
  std::vector<double>  props(20, 0.0);
  props[0] = degree+1;
  props[2] = 0;
  props[3] = 1.0;
  props[4] = 0.01;
  props[5] = 1.0;

  props[8] = 0.0; // SUPG
  props[9] = 1.0; // PSPG
  props[10] = 0.0; // LSIC

  hbs.SetFluidProperties(props);
  
  hbs.SetControl(0, tol1, rho1);


  hbs.setNdof(1);
  
  int parm[5], ii, iter, niter=10, solv=1;
  parm[0] = 2;
  if(solv == 4)
    parm[0] = 2;

  if(solv == 6)
    parm[0] = 4;

 
  hbs.prepareInputData();
  
  cout << " Model preparation successful " << endl;

  hbs.setSolver(solv, parm);

  cout << " Solver preparation successful " << endl;
  
  hbs.setTimeParam();
  
  
  timeFunction.add(new TimeFunction);

  timeFunction[0].fct.add(new TimeFunctionCore);

  timeFunction[0].fct[0].t0 = 0.0;
  timeFunction[0].fct[0].t1 = 1.0;
  timeFunction[0].fct[0].p[0] = 1.0;
  timeFunction[0].fct[0].p[1] = 0.0;
  timeFunction[0].fct[0].p[2] = 0.0;
  timeFunction[0].fct[0].p[3] = 0.0;
  timeFunction[0].fct[0].p[4] = 0.0;
  timeFunction[0].fct[0].p[5] = 0.0;

  timeFunction[0].fct[0].tp = 1.0e+30;


  for(ii=0; ii<timeFunction.n; ii++)
    timeFunction[ii].update();

  hbs.timeUpdate();

  for(iter=0; iter<niter; iter++)
  {
    //cout << " aaaaaaaaaaa " << endl;
    hbs.calcStiffnessAndResidual();

    //cout << " bbbbbbbbbbb " << endl;
    if( hbs.converged() )
      break;

    //cout << " ccccccccccc " << endl;
    hbs.factoriseSolveAndUpdate();

    //cout << " ddddddddddd " << endl;
    hbs.updateIterStep();
  }

  hbs.computeElementErrors(0);
  hbs.computeElementErrors(1);

  //hbs.computeConditionNumber();

  //hbs.printData(8, 0);
  int resln[3];
  resln[0]= resln[1] =resln[2]=1;

  hbs.postProcessFlow(0, 0, 10, 1, 0.1, 1.0, resln);


  printf("\n\n\n\n\n\n\n\n\n\n\n");


  return 0;
}
*/



int main(int argc, char* argv[])
{
  PetscErrorCode ierr;

  //PetscInitialize(NULL,NULL,(char *)0,NULL);

  PetscInitialize(NULL, NULL, "petsc_options.dat", NULL);

  PetscInt  n_mpi_procs, this_mpi_proc;

  MPI_Comm_size(MPI_COMM_WORLD, &n_mpi_procs);

  MPI_Comm_rank(MPI_COMM_WORLD, &this_mpi_proc);


 // Laplace/Stokes problem in 2D

  double  tol1 = 1.0e-6;
  double  rho1 = 0.8;

  int  resln[]={1,1,1};  
  int  ee, ne = 10, degree=1;

  double  PEN  = 1.0e4;

  bool      isNitsche = true;
  double  NitscheFact = 1.0;


  if(argc == 0)
    cerr << " Input data " << endl;
  else
  {
    ne  = atoi(argv[1]);
    degree  = atoi(argv[2]);
    //PEN = atof(argv[3]);
    //isNitsche = (atoi(argv[4]) == 1);
    //NitscheFact = atof(argv[5]);
  }


    HBSplineCutFEM  hbs;

    cout << "Hello, world! " << endl;

    hbs.SetDimension(2);

    hbs.setNdof(1);

    hbs.SetOrigin(0.0, 0.0);

    hbs.SetGridDimensions(1.0, 1.0);
  
    hbs.SetBSplineDegree(degree, degree);
  
    hbs.SetNumberOfElements(ne, ne);


    hbs.AddDirichletBCs(1, 1, 0.0, PEN, isNitsche, NitscheFact);
    hbs.AddDirichletBCs(2, 1, 0.0, PEN, isNitsche, NitscheFact);
    hbs.AddDirichletBCs(3, 1, 0.0, PEN, isNitsche, NitscheFact);
    hbs.AddDirichletBCs(4, 1, 0.0, PEN, isNitsche, NitscheFact);

    //hbs.AddDirichletBCs(1, 2, 0.0, PEN, isNitsche, NitscheFact);
    //hbs.AddDirichletBCs(2, 2, 0.0, PEN, isNitsche, NitscheFact);
    //hbs.AddDirichletBCs(3, 2, 0.0, PEN, isNitsche, NitscheFact);
    //hbs.AddDirichletBCs(4, 2, 0.0, PEN, isNitsche, NitscheFact);

    //hbs.AddPointBCs(0.0, 0.0, 1.0, 3, 0.0, 100.0);
  
    std::vector<double>  props(20, 0.0);

    props[0] = degree+1;
    props[2] = 0;
    props[3] = 1.0;
    props[4] = 1.0;
    props[5] = 1.0;

    props[8] = 0.0; // SUPG
    props[9] = 0.0; // PSPG
    props[10] = 0.0; // LSIC

    hbs.SetFluidProperties(props);

    hbs.SetControl(0, tol1, rho1);

    int parm[5], ii, iter, niter=5, solv=8;

    parm[0] = 1;
    if(solv == 4)
      parm[0] = 1;

    if(solv == 6)
      parm[0] = 4;

 
    hbs.prepareInputData();
  
    cout << " Model preparation successful " << endl;

    hbs.setSolver(solv, parm);

    cout << " Solver preparation successful " << endl;
  
    hbs.setTimeParam();
  
    if(this_mpi_proc == 0)
      hbs.plotGeom(0, 0, 0, 0, resln);

    timeFunction.add(new TimeFunction);

    timeFunction[0].fct.add(new TimeFunctionCore);

    timeFunction[0].fct[0].t0 = 0.0;
    timeFunction[0].fct[0].t1 = 100.0;

    timeFunction[0].fct[0].p[0] = 1.0;
    timeFunction[0].fct[0].p[1] = 0.0;
    timeFunction[0].fct[0].p[2] = 0.0;
    timeFunction[0].fct[0].p[3] = 0.0;
    timeFunction[0].fct[0].p[4] = 0.0;
    timeFunction[0].fct[0].p[5] = 0.0;
    timeFunction[0].fct[0].p[6] = 0.0;
    timeFunction[0].fct[0].p[7] = 0.0;
    timeFunction[0].fct[0].p[8] = 0.0;


    timeFunction[0].fct[0].tp = 1.0e+30;


    for(ii=0; ii<timeFunction.n; ii++)
      timeFunction[ii].update();

    hbs.timeUpdate();

    for(iter=0; iter<niter; iter++)
    {
      //cout << " aaaaaaaaaaa " << endl;
      hbs.calcStiffnessAndResidual();

      //cout << " bbbbbbbbbbb " << endl;
      if( hbs.converged() )
        break;

      //cout << " ccccccccccc " << endl;
      hbs.factoriseSolveAndUpdate();

      //cout << " ddddddddddd " << endl;
      hbs.updateIterStep();
    }

    if(this_mpi_proc == 0)
      hbs.postProcessFlow(0, 0, 1, 0, 0.0, 1.0, resln);

    //hbs.computeElementErrors(0);
    //hbs.computeElementErrors(1);
    //hbs.computeElementErrors(2);
    //hbs.computeElementErrors(3);

    //hbs.computeConditionNumber();

    //hbs.printData(8, 0);

    //printf("\n\n\n\n\n\n\n\n\n\n\n");
  //}

    ierr = PetscFinalize();//CHKERRQ(ierr);

    cout << " \n\n\n Propgram successful ... \n\n\n " << endl;


  return 0;
}















