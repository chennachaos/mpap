

#include "headersBasic.h"

#include "HBSplineFEM.h"
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



/*
int main(int argc, char* argv[])
{
 // Laplace/Stokes problem in 2D
  
  
  double  PEN  ;
  double  tol1 = 1.0e-6;
  double  rho1 = 0.8;
  double  NitscheFact;
  bool  isNitsche;
  
  int ne = 10, degree=1;

  PEN  = 1.0e6;

  isNitsche = false;
  isNitsche = true;

  NitscheFact = 1.0;


  
    if(argc == 0)
      cerr << " Input data " << endl;
    else
    {
       ne  = atoi(argv[1]);
       degree  = atoi(argv[2]);
    }
    
    //PEN /= (1.0/ne);


  HBSplineFEM  hbs;

  cout << "Hello, world! " << endl;
  
  hbs.SetDimension(2);

  hbs.SetOrigin(0.0, 0.0);
  
  hbs.SetGridDimensions(1.0, 1.0);
  
  hbs.SetBSplineDegree(degree, degree);
  
  hbs.SetNumberOfElements(ne, ne);



  hbs.AddDirichletBCs(1, 1, 0.0, PEN, isNitsche, NitscheFact);
  hbs.AddDirichletBCs(2, 1, 0.0, PEN, isNitsche, NitscheFact);
  hbs.AddDirichletBCs(3, 1, 0.0, PEN, isNitsche, NitscheFact);
  hbs.AddDirichletBCs(4, 1, 0.0, PEN, isNitsche, NitscheFact);

  hbs.AddDirichletBCs(1, 2, 0.0, PEN, isNitsche, NitscheFact);
  hbs.AddDirichletBCs(2, 2, 0.0, PEN, isNitsche, NitscheFact);
  hbs.AddDirichletBCs(3, 2, 0.0, PEN, isNitsche, NitscheFact);
  hbs.AddDirichletBCs(4, 2, 0.0, PEN, isNitsche, NitscheFact);


  hbs.AddPointBCs(0.0, 0.0, 1.0, 3, 0.0, 100.0);
  
  std::vector<double>  props(20, 0.0);
  props[0] = degree+1;
  props[2] = 0;
  props[3] = 1.0;
  props[4] = 1.0;

  props[8] = 0.0; // SUPG
  props[9] = 1.0; // PSPG
  props[10] = 0.0; // LSIC

  hbs.SetFluidProperties(props);
  
  hbs.SetControl(0, tol1, rho1);


  hbs.setNdof(3);
  
  int parm[5], ii, iter, niter=10, solv=6;
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
  hbs.computeElementErrors(2);
  hbs.computeElementErrors(3);


  //hbs.computeConditionNumber();

  //hbs.printData(8, 0);

  return 0;
}
*/




/*
int main()
{

  double  PEN  = 1.0e6;
  double  tol1 = 1.0e-6;
  double  rho1 = 0.8;

  int parm[5], ii, iter, niter=10, ne, ee;
  parm[0] = 4;


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
  

  for(ee=0; ee<2; ee++)
  {
    ne = (ee+1)*10;

    //TreeNode<2>::nodecount = 0;

    HBSplineFEM  hbs;

    hbs.SetDimension(2);

    hbs.SetOrigin(0.0, 0.0);
  
    hbs.SetGridDimensions(1.0, 1.0);
  
    hbs.SetBSplineDegree(1, 1);
  
    hbs.SetNumberOfElements(ne, ne);

    hbs.AddDirichletBCs(1, 1, 0.0, PEN);
    hbs.AddDirichletBCs(2, 1, 0.0, PEN);
    hbs.AddDirichletBCs(3, 1, 0.0, PEN);
    hbs.AddDirichletBCs(4, 1, 0.0, PEN);

    std::vector<double>  props(20, 0.0);
    props[0] = 2; // nGP
    props[2] = 0;
    props[3] = 1.0;
    props[4] = 1.0;

    props[8] = 0.0; // SUPG
    props[9] = 1.0; // PSPG
    props[10] = 0.0; // LSIC

    hbs.SetFluidProperties(props);
  
    hbs.SetControl(0, tol1, rho1);

    hbs.setNdof(1);
  
    hbs.prepareInputData();
  
    cout << " Model preparation successful " << endl;

    hbs.setSolver(1, parm);

    cout << " Solver preparation successful " << endl;
  
    hbs.setTimeParam();
  
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
  }

  return 0;
}
*/






int main(int argc, char* argv[])
{
 // Poisson equation 1D


  double  tol1 = 1.0e-6;
  double  rho1 = 0.8;

  
  int  ee, ne = 10, degree=1;

  double  PEN  = 1.0e6;

  //bool  isNitsche = false;
  bool  isNitsche = true;

  double  NitscheFact = 1.0;


  if(argc == 0)
    cerr << " Input data " << endl;
  else
  {
    ne  = atoi(argv[1]);
    degree  = atoi(argv[2]);
    PEN = atof(argv[3]);
    isNitsche = (atoi(argv[4]) == 1);
    NitscheFact = atof(argv[5]);
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
  props[4] = 1.0;
  props[5] = 1.0;
  props[6] = 1.0;


  props[8] = 0.0; // SUPG
  props[9] = 0.0; // PSPG
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

  hbs.computeElementErrors(0);
  hbs.computeElementErrors(2);

  hbs.computeConditionNumber();

  //hbs.printData(8, 0);
  //int resln[3];
  //resln[0]= resln[1] =resln[2]=1;

  //hbs.postProcessFlow(0, 0, 10, 1, 0.1, 1.0, resln);

  printf("\n\n\n\n\n\n\n\n\n\n\n");



  return 0;
}









