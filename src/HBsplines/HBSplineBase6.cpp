
#include "HBSplineBase.h"

#include "SolverEigen.h"
#include "SolverMA41Eigen.h"
#include "SolverPardisoEigen.h"

#include "SolverPetsc.h"
#include "SolverPardisoPetsc.h"
#include "MpapTime.h"

extern  MpapTime  mpapTime;



int  HBSplineBase::fsi_monolithic_fixedpoint_forcePred(int max_iter, double tol_local)
{
  //cout << " HBSplineBase::fsi_monolithic_fixedpoint_forcePred() " << endl;
  PetscPrintf(MPI_COMM_WORLD, " HBSplineBase::fsi_monolithic_fixedpoint_forcePred() \n");

  //
  // Thesis by Sudhakar Yogaraj @TUM, Germany
  // An embedded interface finite element method for fluid-structure-fracture interaction
  // 

  // iterate
  //   1.) predict force on the solid ---> fSP
  //   2.) solve solid problem to get structural displacement ---> dS
  //   3.) Update the interface location
  //   4.) Compute interface velocity
  //   5.) Solve the fluid problem
  //   6.) Compute interface force
  //   7.) IF, converged
  //          update the solution
  //       ELSE
  //          compute relaxation parameter
  //          compute relaxed force
  //       END
  // end iterate

  VectorXd  resiIntf, resiIntfPrev, forceIntf, forceIntfPrev;
  VectorXd  vecTemp;
  double  relaxPara=0.01, relaxParaPrev, normTemp;

  bool  converg=false;
  int  iter=0, bb, resln[]={1,1,1};

  resiIntfPrev = ImmersedBodyObjects[bb]->SolnData.force ;
  //resiIntfPrev.setZero();
  resiIntf = resiIntfPrev;

  while(!converg || (iter < 2))
  {
    //   1.) predict force on the solid ---> fSP

    for(bb=0; bb<nImmSolids; bb++)
    {
      ImmersedBodyObjects[bb]->SolnData.forceCur = ImmersedBodyObjects[bb]->SolnData.force ;

      //ImmersedBodyObjects[bb]->SolnData.forceCur = 2.0*ImmersedBodyObjects[bb]->SolnData.force - ImmersedBodyObjects[bb]->SolnData.forcePrev;
    }

    //   2.) solve solid problem to get structural displacement ---> dS
    //   3.) Update the interface location
    //   4.) Compute interface velocity

    solveSolidProblem();

    updateImmersedPointPositions();

    updateIterStep();

    //   5.) Solve the fluid problem
    solveFluidProblem();

    //   6.) Compute interface force
    //   7.) IF, converged
    //          update the solution
    //       ELSE
    //          compute relaxation parameter
    //          compute relaxed force
    //       END

    postProcessFlow(1, 1, 10, 1, 0.0, 1.0, resln);

    for(bb=0; bb<nImmSolids; bb++)
    {
      computeTotalForce(bb);

      ImmersedBodyObjects[bb]->updateForce(&(totalForce(0)));

      resiIntfPrev = resiIntf;

      resiIntf = ImmersedBodyObjects[bb]->SolnData.force - ImmersedBodyObjects[bb]->SolnData.forcePrev;

      normTemp = resiIntf.norm()/sqrt(resiIntf.rows());

      PetscPrintf(MPI_COMM_WORLD, " Force convergence ...  Iter = %5d \t \t %11.4E\n", iter, normTemp);

      if(normTemp <= tol_local)
      {
        converg = true;
      }
      else
      {
        converg = false;

        vecTemp = resiIntf - resiIntfPrev;

        normTemp = vecTemp.squaredNorm();

        relaxParaPrev = relaxPara;
        
        if(iter > 1)
        {
          //ImmersedBodyObjects[bb]->perform_Aitken_accelerator_force();
          //relaxPara = min(-relaxPara*(resiIntfPrev.dot(vecTemp)/normTemp), relaxParaPrev);
          relaxPara = -relaxPara*(resiIntfPrev.dot(vecTemp)/normTemp);
        }
        //else
        //{
          //relaxPara = 0.01;
        //}

        cout << " iter = " << iter << " \t relaxPara = " << relaxPara << endl;

        ImmersedBodyObjects[bb]->SolnData.force  = relaxPara*ImmersedBodyObjects[bb]->SolnData.force;
        ImmersedBodyObjects[bb]->SolnData.force  += (1.0-relaxPara)*ImmersedBodyObjects[bb]->SolnData.forcePrev;
      }

      ImmersedBodyObjects[bb]->SolnData.forcePrev3 = ImmersedBodyObjects[bb]->SolnData.forcePrev2;
      ImmersedBodyObjects[bb]->SolnData.forcePrev2 = ImmersedBodyObjects[bb]->SolnData.forcePrev;
      ImmersedBodyObjects[bb]->SolnData.forcePrev  = ImmersedBodyObjects[bb]->SolnData.force;

    } // for(bb=0; bb<nImmSolids; bb++)

    iter++;
  }
  
  return 1;
}

/*

    for(bb=0; bb<nImmSolids; bb++)
    {
      computeTotalForce(bb);

      ImmersedBodyObjects[bb]->updateForce(&(totalForce(0)));

      resiIntfPrev = resiIntf;

      resiIntf = ImmersedBodyObjects[bb]->SolnData.force - ImmersedBodyObjects[bb]->SolnData.forcePrev;

      vecTemp = resiIntf - resiIntfPrev;

      normTemp = vecTemp.norm();

      if(normTemp <= tolTemp)
      {
        converg = true;
      }
      else
      {
        converg = false;

        normTemp = normTemp*normTemp;

        w  = -wP*resiIntfPrev.dot(vecTemp)/normTemp;
        wP = w;

        cout << " w = " << w << endl;

        ImmersedBodyObjects[bb]->SolnData.force = w*ImmersedBodyObjects[bb]->SolnData.force + (1.0-w)*ImmersedBodyObjects[bb]->SolnData.forcePrev;

        ImmersedBodyObjects[bb]->SolnData.forceCur = ImmersedBodyObjects[bb]->SolnData.force;
      }
    } // for(bb=0; bb<nImmSolids; bb++)
*/




int  HBSplineBase::fsi_monolithic_fixedpoint_dispPred(int max_iter, double tol_local)
{
  //cout << " HBSplineBase::fsi_staggered_displacement_predictor() " << endl;

  //
  // Ulrich Küttler and Wolfgang A. Wall
  // Fixed-point fluid–structure interaction solvers with dynamic relaxation
  // Comput Mech (2008) 43:61–72
  // DOI 10.1007/s00466-008-0255-5
  //
  // iterate
  //   Step 1: predict interface displacement ---> dP
  //   Step 2: solve the mesh problem  ---> vI
  //   Step 3: solve fluid problem with the interface velocity obtained from Step 2 ---> fS
  //   Step 4: a.)solve solid problem with the force obtained from Step 3 ---> dS
  //           b.) check for convergence.
  //           c.) update the displacement using Aitken accelerator
  // end

  // or

  // A fixed-grid b-spline finite element technique for fluid-structure interaction
  // Ruberg and Cirak, IJNME, 74:623:660, 2013.

  //   1.) predict the displacement and velocity of the solid/interface ---> dS and vS
  // iterate
  //   2.) Update the interface location
  //   3.) Solve the fluid problem
  //   4.) Compute the interface force
  //   5.) solve solid problem to get structural displacement ---> dS
  //   6.) Check convergence. IF not converged THEN go to Step 2.


  VectorXd  resiIntf, resiIntfPrev, forceIntf, forceIntfPrev;
  VectorXd  vecTemp;
  double  relaxPara=0.1, relaxParaPrev, normTemp;

  bool  converg=false;
  int  iter=0, bb, resln[]={1,1,1};

  resiIntfPrev = ImmersedBodyObjects[bb]->SolnData.var1 ;
  resiIntfPrev.setZero();
  resiIntf = resiIntfPrev;

  //   1.) predict the displacement and velocity of the solid/interface ---> dS and vS

    for(bb=0; bb<nImmSolids; bb++)
    {
      ImmersedBodyObjects[bb]->initialise_solid_state();
    }


  //while( (iter < max_iter) || converg)
  while( !converg )
  {
    //   2.) Update the interface location
    //   3.) Solve the fluid problem

    updateImmersedPointPositions();

    updateIterStep();

    solveFluidProblem();

    //postProcessFlow(1, 1, 10, 1, 0.0, 1.0, resln);

    //   4.) Compute the interface force

    for(bb=0; bb<nImmSolids; bb++)
    {
      computeTotalForce(bb);

      ImmersedBodyObjects[bb]->updateForce(&(totalForce(0)));

    } // for(bb=0; bb<nImmSolids; bb++)


    //   5.) solve solid problem to get structural displacement ---> dS

    solveSolidProblem();

    //   6.) Check convergence. IF not converged THEN go to Step 2.

    for(bb=0; bb<nImmSolids; bb++)
    {
      resiIntfPrev = resiIntf;

      resiIntf = ImmersedBodyObjects[bb]->SolnData.var1 - ImmersedBodyObjects[bb]->SolnData.var1PrevIter;

      normTemp = resiIntf.norm()/sqrt(resiIntf.rows());

      PetscPrintf(MPI_COMM_WORLD, " Displacement convergence ...  Iter = %5d \t \t %11.4E\n", iter, normTemp);

      if(normTemp <= tol_local)
      {
        converg = true;
      }
      else
      {
        converg = false;

        vecTemp = resiIntf - resiIntfPrev;

        normTemp = vecTemp.squaredNorm();

        relaxParaPrev = relaxPara;

        if(iter > 1)
        {
          //ImmersedBodyObjects[bb]->perform_Aitken_accelerator_force();
          //relaxPara = min(-relaxPara*(resiIntfPrev.dot(vecTemp)/normTemp), relaxParaPrev);
          relaxPara = -relaxPara*(resiIntfPrev.dot(vecTemp)/normTemp);
        }
        //else
        //{
          //relaxPara = 0.01;
        //}

        cout << " iter = " << iter << " \t relaxPara = " << relaxPara << endl;

        ImmersedBodyObjects[bb]->SolnData.var1  = relaxPara*ImmersedBodyObjects[bb]->SolnData.var1;
        ImmersedBodyObjects[bb]->SolnData.var1  += (1.0-relaxPara)*ImmersedBodyObjects[bb]->SolnData.var1PrevIter;
      }

      ImmersedBodyObjects[bb]->SolnData.var1PrevIter = ImmersedBodyObjects[bb]->SolnData.var1;
      ImmersedBodyObjects[bb]->updateIterStep();
    } // for(bb=0; bb<nImmSolids; bb++)

    iter++;
  }


    char        tmp[100];
    MyString    tmpStr;

    sprintf(tmp," \t %6d", iter);

    tmpStr.append(tmp);
    prgWriteToTFile(tmpStr);

  return 1;
}





int  HBSplineBase::fsi_staggered_force_predictor(int max_iter, double tol_local)
{
  PetscPrintf(MPI_COMM_WORLD, "\n ---------------------------------------------------- \n");
  PetscPrintf(MPI_COMM_WORLD, "\n\n HBSplineBase::fsi_staggered_force_predictor()      \n");
  PetscPrintf(MPI_COMM_WORLD, "\n ---------------------------------------------------- \n\n\n");

  //
  // Staggered scheme by Dettmer and Peric
  // A new staggered scheme for fluid-structure interaction
  // 

  //   1.) predict force on the solid ---> fSP
  //   2.) solve solid problem to get structural displacement ---> dS
  //   3.) Update the interface location
  //   4.) Compute interface velocity
  //   5.) Solve the fluid problem
  //   6.) Compute the interface force
  //   7.) Relax interface force


  VectorXd  resiIntf, resiIntfPrev, forceIntf, forceIntfPrev;
  VectorXd  vecTemp;
  double  normTemp, tolTemp=1.0e-4;

  bool  converg=false;
  int  iter=0, bb, resln[]={1,1,1};

  int  predType = (int) ImmersedBodyObjects[0]->SolnData.stagParams[1];

  double  q1, q2, q3, q4, fact;
  double  knp1 = mpapTime.dt/mpapTime.dtPrev;

  //double  z = (1.0+beta+beta*rho-2.0*sqrt(beta+beta*rho))/(1.0-beta);
  //double  z=0.5;

  switch(predType)
  {
      case 1:
        q1 =  1.0;  q2 = 0.0;  q3 =  0.0;  q4 =  0.0;
        break;

      case 2:
        //q1 =  2.0;  q2 = -1.0;  q3 =  0.0;  q4 =  0.0;
        q1 = 1.0+knp1; q2 = -knp1;  q3 =  0.0;  q4 =  0.0;
        break;

      case 3:
        q1 =  3.0;  q2 = -3.0;  q3 =  1.0;  q4 =  0.0;
        //q1 =  2.0+z;  q2 = -1.0-2.0*z;  q3 =  z;  q4 =  0.0;
        //q1 =  2.5;  q2 = -2.0;  q3 =  0.5;  q4 =  0.0;
        break;

      case 4:
        q1 =  4.0;  q2 = -6.0;  q3 =  4.0;  q4 = -1.0;
        //q1 =  104.0/35.0;  q2 = -114.0/35.0;  q3 =  56.0/35.0;  q4 = -11.0/35.0;
        break;

      default:
        cout << " unknown value for 'predType' " << endl;
        break;
  }

  double beta = ImmersedBodyObjects[0]->SolnData.stagParams[2];

  PetscPrintf(MPI_COMM_WORLD, " Staggered scheme parameters: ... \n\n");
  PetscPrintf(MPI_COMM_WORLD, " Predictor type              = %d \n", predType);
  PetscPrintf(MPI_COMM_WORLD, " dt/dtPrev                   = %f \n", knp1);
  PetscPrintf(MPI_COMM_WORLD, " Relaxation parameter (beta) = %f \n\n", beta);


  PetscPrintf(MPI_COMM_WORLD, "\n timeUpdate   \n");
  PetscPrintf(MPI_COMM_WORLD, " -------------  \n");

  timeUpdate();


  //   1.) predict force on the solid ---> fSP
  for(bb=0; bb<nImmSolids; bb++)
  {
    ImmersedBodyObjects[bb]->SolnData.force = q1*ImmersedBodyObjects[bb]->SolnData.forcePrev + q2*ImmersedBodyObjects[bb]->SolnData.forcePrev2;
  }


  for(iter=0; iter<max_iter; iter++)
  {
    PetscPrintf(MPI_COMM_WORLD, " Global iteration = %d \n", iter);
    PetscPrintf(MPI_COMM_WORLD, " ------------------------- \n");

    for(bb=0; bb<nImmSolids; bb++)
    {
      double  af =  ImmersedBodyObjects[bb]->SolnData.td[2];

      ImmersedBodyObjects[bb]->SolnData.forceCur = af*ImmersedBodyObjects[bb]->SolnData.force + (1.0-af)*ImmersedBodyObjects[bb]->SolnData.forcePrev;

      ImmersedBodyObjects[bb]->SolnData.var1PrevIteration = ImmersedBodyObjects[bb]->SolnData.var1;
    }

    //   2.) solve solid problem to get structural displacement ---> dS

    PetscPrintf(MPI_COMM_WORLD, "\n solveSolidProblem    \n");
    PetscPrintf(MPI_COMM_WORLD, " ---------------------  \n");

    solveSolidProblem();

    for(bb=0; bb<nImmSolids; bb++)
    {
      VectorXd  vecTemp = ImmersedBodyObjects[bb]->SolnData.var1PrevIteration - ImmersedBodyObjects[bb]->SolnData.var1;
      normTemp = vecTemp.norm();

      PetscPrintf(MPI_COMM_WORLD, " Displacement convergence ...  Iter = %5d \t \t %11.4E\n", iter, normTemp);
    }

    //   3.) Update the interface location
    //   4.) Compute interface velocity

    PetscPrintf(MPI_COMM_WORLD, "\n updateImmersedPointPositions  \n");
    PetscPrintf(MPI_COMM_WORLD, " ------------------------------  \n");
    updateImmersedPointPositions();

    PetscPrintf(MPI_COMM_WORLD, "\n updateIterStep       \n");
    PetscPrintf(MPI_COMM_WORLD, " ---------------------  \n");
    updateIterStep();

    //   5.) Solve the fluid problem
    PetscPrintf(MPI_COMM_WORLD, "\n solveFluidProblem    \n");
    PetscPrintf(MPI_COMM_WORLD, " ---------------------  \n");
    solveFluidProblem();

    //   6.) Compute the interface force
    //   7.) Relax interface force

    //postProcessFlow(1, 1, 10, 1, 0.0, 1.0, resln);

    PetscPrintf(MPI_COMM_WORLD, "\n relax force \n");
    PetscPrintf(MPI_COMM_WORLD, " ---------------------  \n");
    for(bb=0; bb<nImmSolids; bb++)
    {
      //cout << " computeTotalForce() " << endl;
      // This is for CutFEM
      computeTotalForce(bb);

      //PetscPrintf(MPI_COMM_WORLD, " update force \n");
      // This is for CutFEM
      ImmersedBodyObjects[bb]->updateForce(&(totalForce(0)));
      // This is for FDM
      //ImmersedBodyObjects[bb]->updateForce();

      ImmersedBodyObjects[bb]->SolnData.forcePrevIteration  =  ImmersedBodyObjects[bb]->SolnData.force;

      ImmersedBodyObjects[bb]->SolnData.force  =  beta*ImmersedBodyObjects[bb]->SolnData.forceTemp + (1.0-beta)*ImmersedBodyObjects[bb]->SolnData.force;

      VectorXd  vecTemp = ImmersedBodyObjects[bb]->SolnData.force - ImmersedBodyObjects[bb]->SolnData.forcePrevIteration; 
      normTemp = vecTemp.norm();

      PetscPrintf(MPI_COMM_WORLD, " Force convergence ...  Iter = %5d \t \t %11.4E\n", iter, normTemp);
    } // for(bb=0; bb<nImmSolids; bb++)
    //PetscPrintf(MPI_COMM_WORLD, " relax force DONE \n");
  }

    //for(bb=0; bb<nImmSolids; bb++)
    //{
      //ImmersedBodyObjects[bb]->SolnData.forcePrev3 = ImmersedBodyObjects[bb]->SolnData.forcePrev2;
      //ImmersedBodyObjects[bb]->SolnData.forcePrev2 = ImmersedBodyObjects[bb]->SolnData.forcePrev;
      //ImmersedBodyObjects[bb]->SolnData.forcePrev  = ImmersedBodyObjects[bb]->SolnData.force;
    //}

    PetscPrintf(MPI_COMM_WORLD, "\n ---------------------------------------------------- \n");
    PetscPrintf(MPI_COMM_WORLD, "\n ---------------------------------------------------- \n");
    PetscPrintf(MPI_COMM_WORLD, "\n ---------------------------------------------------- \n\n");

    return 0;
}


int  HBSplineBase::fsi_staggered_displacement_predictor(int max_iter, double tol_local)
{
  cout << " HBSplineBase::fsi_staggered_displacement_predictor() " << endl;

  //
  // 
  //   1.) predict the displacement of the solid ---> dSP
  //   2.) Update the interface location
  //   3.) Compute interface velocity
  //   4.) Solve the fluid problem
  //   5.) Compute the interface force
  //   6.) solve solid problem to get structural displacement ---> dS
  //   7.) Relax interface displacement


  VectorXd  resiIntf, resiIntfPrev, forceIntf, forceIntfPrev;
  VectorXd  vecTemp;
  double  relaxPara=0.5, relaxParaPrev, normTemp, tolTemp=1.0e-4;

  bool  converg=false;
  int  iter=0, bb, resln[]={1,1,1};

  resiIntfPrev = ImmersedBodyObjects[bb]->SolnData.force ;
  resiIntfPrev.setZero();
  resiIntf = resiIntfPrev;

  while(iter < 1)
  {
    //   1.) predict the displacement of the solid ---> dSP

    for(bb=0; bb<nImmSolids; bb++)
    {
      //ImmersedBodyObjects[bb]->SolnData.forceCur = ImmersedBodyObjects[bb]->SolnData.force ;

      ImmersedBodyObjects[bb]->SolnData.forceCur = 2.0*ImmersedBodyObjects[bb]->SolnData.force - ImmersedBodyObjects[bb]->SolnData.forcePrev;
    }

    //   2.) solve solid problem to get structural displacement ---> dS
    //   3.) Update the interface location
    //   4.) Compute interface velocity

    solveSolidProblem();

    updateImmersedPointPositions();

    updateIterStep();

    //   5.) Solve the fluid problem
    solveFluidProblem();

    //   6.) Compute the interface force
    //   7.) Relax interface force

    postProcessFlow(1, 1, 10, 1, 0.0, 1.0, resln);

    for(bb=0; bb<nImmSolids; bb++)
    {
      computeTotalForce(bb);

      ImmersedBodyObjects[bb]->updateForce(&(totalForce(0)));

      //PetscPrintf(MPI_COMM_WORLD, " Force convergence ...  Iter = %5d \t \t %11.4E\n", iter, normTemp);

      //ImmersedBodyObjects[bb]->SolnData.forcePrev3 = ImmersedBodyObjects[bb]->SolnData.forcePrev2;
      //ImmersedBodyObjects[bb]->SolnData.forcePrev2 = ImmersedBodyObjects[bb]->SolnData.forcePrev;
      //ImmersedBodyObjects[bb]->SolnData.forcePrev  = ImmersedBodyObjects[bb]->SolnData.force;
    } // for(bb=0; bb<nImmSolids; bb++)

    iter++;
  }

  return 1;
}



int HBSplineBase::solveFluidProblem()
{
  //printf("\n Solving HBSplineCutFEM::solveFluidProblem() \n");
  //printf("\n External force norm = %12.6E \n", forceCur.norm());

  for(int ii=0;ii<10;ii++)
  {
    calcStiffnessAndResidual(1, 1, 1);

    PetscPrintf(MPI_COMM_WORLD, "\t %5d \t %12.6E\n", ii, rNorm);

    if(converged())
      break;

    //cout << " Fluid ... factoriseSolveAndUpdate " << endl;
    factoriseSolveAndUpdate();
    //cout << " Fluid ... factoriseSolveAndUpdate " << endl;

    //cout << " Fluid ... updateIterStep " << endl;
    updateIterStep();
    //cout << " Fluid ... updateIterStep " << endl;
  }

  PetscPrintf(MPI_COMM_WORLD, "\n Solving HBSplineCutFEM::solveFluidProblem() ..... DONE \n\n");

  return 0;
}

