

#include "HBSplineCutFEM.h"

#include "DataBlockTemplate.h"
#include "SolverEigen.h"
#include "ComputerTime.h"
#include "MpapTime.h"
#include "NurbsUtilities.h"

#include "BasisFunctionsLagrange.h"
#include "myDataIntegrateCutFEM.h"
#include "ImmersedIntegrationElement.h"
#include "QuadratureUtil.h"
#include "typedefs.h"

#include <cmath>
#include <omp.h>


extern ComputerTime       computerTime;
extern MpapTime mpapTime;



void HBSplineCutFEM::setInitialConditions()
{
    double* tmp;

    solverEigen->zeroMtx();

    for(int ee=0;ee<elems.size();ee++)
    {
       if( !(elems[ee]->IsGhost()) &&  elems[ee]->IsLeaf() )
       {
           //cout << " elems[ee]->GetID() " <<  elems[ee]->GetID() << '\t' <<  elems[ee]->GetLevel() << endl;

           elems[ee]->resetMatrixAndVector();
           elems[ee]->setInitialProfile();
           //elems[ee]->AssembleMatrixAndVector(1, solver->mtx, &(rhsVec(0)));
       }
    }

    solverEigen->currentStatus = ASSEMBLY_OK;

    factoriseSolveAndUpdate();

    //rhsVec = rhsVec * (1.0/rhsVec.maxCoeff());
    //

    //
    int resln1[3]; resln1[0] = resln1[1] = 5;

    //postProcess2D(1, 1, 10, 1, 0.0, 1.0, resln1);
    //postProcess1D(1, 1, 10, 1, 0.0, 1.0, resln1);
    //

    solnInit = soln;
    SolnData.var1 = solnInit;

   return;
}




void  HBSplineCutFEM::applyBoundaryConditions()
{
    //cout << "     HBSplineCutFEM: applyBoundaryConditions ... STARTED  \n\n";
    //firstIter = true;

    //cout << "     HBSplineCutFEM: applyBoundaryConditions ...DONE  \n\n";
    return;
}



void HBSplineCutFEM::addExternalForces()
{
//  cout << "HBSplineCutFEM::addExternalForces() " << endl;
  // if(firstIter)     calcAndAssyLoadVector(1.0, 0.0);

  return;
}




void HBSplineCutFEM::computeElementErrors(int index)
{
    if(index==5)
    {
      solverEigen->computeConditionNumber();
      return;
    }

    int  ii, ee, count=0, dd, domTemp;
    node  *nd;

    totalError = 0.0;
    
    if(index < 4) // L2 or H1 norm based errors
    {
       for(ee=0; ee<activeElements.size(); ee++)
       {
          nd = elems[activeElements[ee]];

          //for(dd=0; dd<nd->domainNums.size(); dd++)
          //{
            //domTemp = nd->domainNums[dd] ;
            domTemp = nd->GetDomainNumber() ;
            //if(domTemp == 1)
            //{
              nd->calcError(index, domTemp);
            
              //cout << ee << '\t' << elems[ee]->GetError() << endl;
              //totalError += ( elems[ee]->GetError() * elems[ee]->GetError() );
              totalError +=  nd->GetError();
            //}
          //}
       }
       totalError = sqrt(totalError);
    }
    else // gradient based error
    {
       for(ii=0;ii<activeElements.size();ii++)
       {
          nd = elems[activeElements[ii]];

          nd->calcError(index);
            
          totalError += nd->GetError();

          count++;
       }
       totalError /= count;
    }
    
    //printf(" \n\n \t totalError = %12.6E \n\n " , totalError);
    
    if(index < 3)
      printf(" \n\n \t L2 Error = %12.6E \n\n " , totalError);
    else
      printf(" \n\n \t H1 Error = %12.6E \n\n " , totalError);

    char        tmp[50];
    MyString    tmpStr, fname;

    sprintf(tmp,"\n %5d \t %5d \t %12.6E ", CURRENT_LEVEL, totalDOF, totalError);
    tmpStr.append(tmp);

    prgWriteToTFile(tmpStr);
    


    return;
}





void HBSplineCutFEM::computeConditionNumber()
{
  solverEigen->computeConditionNumber();

  return;
}



int HBSplineCutFEM::calcStiffnessAndResidual(int solver_type, bool zeroMtx, bool zeroRes)
{
    //cout << "     HBSplineCutFEM: generating coefficient Matrices ...\n\n";

    //assert( SOLVER_TYPE == SOLVER_TYPE_PETSC );

    char fct[] = "HBSplineCutFEM::calcStiffnessAndResidual";
    computerTime.go(fct);

    int bb, ee, ii, jj, kk, ll, nr, nc, aa, dd, r, c, gp, start=0;
    int pp, dof, dir, domTemp=0;
    double  PENALTY;

    time_t tstart, tend;

    IterNum   = (iterCount == 1);
    
    //cout << " IterNum " << IterNum << endl;

    if(firstIter)
    {
      //solver->ResetPrecondFlag();
      rNorm = -1.0;
    }

    ////////////////////////////////////////////////////////
    //
    // stiffness and residual for the background fluid grid
    //
    ////////////////////////////////////////////////////////

    //cout << " Eigen::getNbThreads() = " << Eigen::nbThreads() << endl;

    tstart = time(0);

    if(SOLVER_TYPE == SOLVER_TYPE_PETSC)
      solverPetsc->zeroMtx();
    else
      solverEigen->zeroMtx();

      //cout << " wwwwwwwwwww " << endl;

      //printVector(&(SolnData.var1(0)), 100);

      /*
      for(ii=0;ii<activeElements.size();ii++)
      {
        node *nd = elems[activeElements[ii]];

        //cout << " nd->GetID() " <<  nd->GetID() << '\t' <<  nd->GetLevel() << endl;

        domTemp = nd->GetDomainNumber() ;

        if( domTemp <= 0 )
        {
          nr = nd->forAssyVec.size();

          //cout << nr << '\t' << nc << endl;
          MatrixXd  Klocal;
          VectorXd  Flocal;

          Klocal = MatrixXd::Zero(nr, nr);
          Flocal = VectorXd::Zero(nr);

          //cout << dd << '\t' << domTemp << endl;

           //cout << " AAAAAAAAAAAAAAAAA " << endl;
          nd->calcStiffnessAndResidualCutFEMFluid(Klocal, Flocal, domTemp);
          //cout << " BBBBBBBBBBBBBBBBB " << endl;

          nd->applyDirichletBCsCutFEMFluid(Klocal, Flocal, domTemp);
          //cout << " DDDDDDDDDDDDDDDDD " << endl;
          nd->applyNeumannBCsCutFEMFluid(Klocal, Flocal, domTemp);
          //cout << " BBBBBBBBBBBBBBBBB " << endl;
          solverPetsc->AssembleMatrixAndVectorCutFEM(start, start, nd->forAssyVec, forAssyCutFEM, Klocal, Flocal);
          //cout << " CCCCCCCCCCCCCCCC " << endl;
        }
      }
      */

      //
      for(ii=0;ii<activeElements.size();ii++)
      {
        node *nd = elems[activeElements[ii]];

        if( nd->get_subdomain_id() == this_mpi_proc )
        {
          domTemp = nd->GetDomainNumber() ;

          //cout << " nd->GetID() " <<  nd->GetID() << '\t' <<  nd->GetLevel() << '\t' << domTemp << endl;

          if( domTemp <= 0 )
          {
            nr = nd->forAssyVec.size();

            //cout << nr << '\t' << nc << endl;
            MatrixXd  Klocal;
            VectorXd  Flocal;

            Klocal = MatrixXd::Zero(nr, nr);
            Flocal = VectorXd::Zero(nr);

            //cout << " AAAAAAAAAAAAAAAAA " << endl;
            nd->calcStiffnessAndResidualGFEM(Klocal, Flocal, domTemp);
            //cout << " BBBBBBBBBBBBBBBBB " << endl;
            nd->applyDirichletBCsGFEM(Klocal, Flocal, domTemp);
            //cout << " DDDDDDDDDDDDDDDDD " << endl;

            //printMatrix(Klocal);
            //printf("\n\n\n");
            //printVector(Flocal);

            solverPetsc->AssembleMatrixAndVectorCutFEM(start, start, nd->forAssyVec, forAssyCutFEM, Klocal, Flocal);
            //cout << " CCCCCCCCCCCCCCCC " << endl;
          }
        }
      }
      //

    //cout << " AAAAAAAAAAAAAAAAA " << endl;
    //printVector(pointBCs[0]);

    myDataIntegrateCutFEM  myData;

    //
    if(pointBCs.size() > 0)
    {
    
      for(pp=0;pp<pointBCs.size();pp++)
      {
        geom[0] = pointBCs[pp][0];
        geom[1] = pointBCs[pp][1];
        geom[2] = pointBCs[pp][2];

        dof = int (pointBCs[pp][3] - 1);

        myData.dir = dof;
        myData.specVal[0] = pointBCs[pp][4];

        PENALTY  = pointBCs[pp][5];
        //printVector(pointBCs[0]);

        myData.PENALTY = PENALTY;
        myData.dvol = 1.0*PENALTY;

        ee = findCellNumber(geom);

        //cout << " ee " << ee << endl;

        geometryToParametric(geom, param);
        
        myData.param = param;
        myData.geom  = geom;

        node *nd;
        nd = elems[ee];

        int nr = nd->forAssyVec.size();// + nd->forAssyVec2.size();

        //cout << nr << '\t' << nc << '\t' << '\t' << PENALTY << endl;

        MatrixXd  Klocal;
        VectorXd  Flocal;

        //Klocal = MatrixXd::Zero(nr, nc);
        //Flocal = VectorXd::Zero(nr);

        myData.K1 = MatrixXd::Zero(nr, nr);
        myData.F1 = VectorXd::Zero(nr);

        //dof, param, spec_val, PENALTY, Klocal, Flocal
        nd->applyBoundaryConditionsAtApoint(myData);
        solverPetsc->AssembleMatrixAndVector(0, 0, nd->forAssyVec, myData.K1, myData.F1);
      }
    }
    //

    //cout << solver->mtx << endl;

    //cout << " rhsVec " << endl;        printVector(&(solver->rhsVec(0)), totalDOF);
    
    //printf("\n rhsVec norm = %12.6E \n", solverEigen->rhsVec.norm());
    
    if(DIM == 2)
      applyInterfaceTerms2D();

    if(DIM == 3)
      applyInterfaceTerms3D();

    //printf("\n rhsVec norm = %12.6E \n", solverEigen->rhsVec.norm());

    //cout << " rhsVec " << endl;

    //printData(1, 1);
    //printData(3, 1);
    //cout << " rhsVec " << endl;        printVector(&(solver->rhsVec(0)), totalDOF);
    //printf("\n rhsVec norm = %12.6E \n", solver->rhsVec.norm());

    firstIter = false;
    rNormPrev = rNorm;

    if(SOLVER_TYPE == SOLVER_TYPE_PETSC)
    {
      VecAssemblyBegin(solverPetsc->rhsVec);
      VecAssemblyEnd(solverPetsc->rhsVec);

      VecNorm(solverPetsc->rhsVec, NORM_2, &rNorm);
      solverPetsc->currentStatus = ASSEMBLY_OK;
    }
    else
    {
      rNorm = solverEigen->rhsVec.norm();
      solverEigen->currentStatus = ASSEMBLY_OK;
    }

    //IterNum   = (iterCount == 1);

    SolnData.firstIter = firstIter;

    if( std::isnan(rNorm) )
    {
      cerr << "  NAN found "  << endl;
      exit(0);
    }

    PetscPrintf(MPI_COMM_WORLD, "  %5d \t %11.4e\n", iterCount, rNorm);

    //COUT << domain.name(this); printf("  %11.4e\n",rNorm);

    ctimCalcStiffRes += computerTime.stop(fct);
    //computerTime.stopAndPrint(fct);
   
    //cout << "     HBSplineCutFEM: generating coefficient Matrices ...DONE  \n\n";

    tend = time(0);
    //printf("HBSplineCutFEM::calcStiffnessAndResidual() took %8.4f second(s) \n ", difftime(tend, tstart) );

    iterCount++;

    return 0;
}






int HBSplineCutFEM::factoriseSolveAndUpdate()
{
  //cout << "     HBSplineCutFEM: solving the matrix system ...  \n\n";
  char fct[] = "HBSplineCutFEM::factoriseSolveAndUpdate";
  computerTime.go(fct);

  time_t tstart, tend;

  int dd, kk, ii;

  //solver->computeConditionNumber();

  //cout << " rhsVec " << endl;        printVector(&(solver->rhsVec(0)), totalDOF);

  tstart = time(0);

  if(SOLVER_TYPE == SOLVER_TYPE_PETSC)
  {
    //VecView(solverPetsc->rhsVec, PETSC_VIEWER_STDOUT_WORLD);

    //cout << " PetscSolver " << endl;

    //tstart = time(0);

    solverPetsc->factoriseAndSolve();

    //tend = time(0);
    //printf("HBSplineCutFEM::factoriseSolveAndUpdate() took %8.4f second(s) \n ", difftime(tend, tstart) );

    //VecView(solver2->soln, PETSC_VIEWER_STDOUT_WORLD);

    /////////////////////////////////////////////////////////////////////////////
    // get the solution vector onto all the processors
    /////////////////////////////////////////////////////////////////////////////
  
    Vec            vec_SEQ;
    VecScatter     ctx;
    PetscScalar *arrayTemp;

    VecScatterCreateToAll(solverPetsc->soln, &ctx, &vec_SEQ);
    VecScatterBegin(ctx, solverPetsc->soln, vec_SEQ, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(ctx, solverPetsc->soln, vec_SEQ, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterDestroy(&ctx);

    VecGetArray(vec_SEQ, &arrayTemp);

    // update solution vector

    int ii, jj, dd, kk=0;
    for(ii=0; ii<nNode; ii++)
    {
      //cout << ii << '\t' << node_map_old_to_new[ii] << '\t' << node_map_new_to_old[ii] << endl;
      //jj = node_map_old_to_new[ii]*ndof;
      jj = ii*ndof;
      for(dd=0; dd<ndof; dd++)
        soln[kk++] = arrayTemp[jj+dd];
    }

    //cout << " fluidDOF = " << fluidDOF << endl;
    for(ii=0; ii<fluidDOF; ii++)
      SolnData.var1[ii] += soln[ii];

    //for(ii=0; ii<fluidDOF; ii++)
      //SolnData.var1[forAssyCutFEM2[ii]] += soln[ii];

    //for(ii=0; ii<fluidDOF; ii++)
      //SolnData.var1[forAssyCutFEM2[node_map_new_to_old[ii]]] += soln[ii];

    VecRestoreArray(vec_SEQ, &arrayTemp);
  }
  else
  {
    solverEigen->factoriseAndSolve();

    tend = time(0);
    //printf("HBSplineCutFEM::factoriseSolveAndUpdate() took %8.4f second(s) \n ", difftime(tend, tstart) );

    //cout << " result " << endl;        printVector(&(solver->soln(0)), totalDOF);

    double  *sln = &(solverEigen->soln[0]);

    for(ii=0; ii<fluidDOF; ii++)
      SolnData.var1[forAssyCutFEM2[ii]] += sln[ii];
  
    if(!STAGGERED)
    {
      cout << " need to update this for monolithic scheme " << endl;
      for(int bb=0; bb<ImmersedBodyObjects.size(); bb++)
      {
        //if(ImmersedBodyObjects[bb]->GetNdof() > 0)
        //{
          //cout << " AAAAAAAAAAA " << bb << endl;
          ImmersedBodyObjects[bb]->updateDisplacement(&(sln[fluidDOF]));
          //cout << " AAAAAAAAAAA " << bb << endl;
        //}
      }
    }
  }

  //printVector(SolnData.var1);
  //printf("\n\n\n");

  //computerTime.stopAndPrint(fct);

  ctimFactSolvUpdt += computerTime.stop(fct);

  //cout << " result " << endl;        printVector(&(soln(0)), totalDOF);

  //cout << "     HBSplineCutFEM: solving the matrix system ...DONE  \n\n";

  //computerTime.stopAndPrint(fct);

  return 0;
}



