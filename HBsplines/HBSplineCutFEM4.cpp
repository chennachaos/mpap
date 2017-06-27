
#include "HBSplineCutFEM.h"
#include "SolverEigen.h"
#include "ComputerTime.h"
#include "MpapTime.h"
#include "myDataIntegrateCutFEM.h"
#include "ImmersedIntegrationElement.h"

#include <chrono>
typedef std::chrono::high_resolution_clock Clock;

extern ComputerTime       computerTime;
extern MpapTime mpapTime;



void HBSplineCutFEM::setInitialConditions()
{
    double* tmp;

    solverEigen->zeroMtx();

    for(int ee=0;ee<elems.size();ee++)
    {
       if( !(elems[ee]->isGhost()) &&  elems[ee]->isLeaf() )
       {
           //cout << " elems[ee]->getID() " <<  elems[ee]->getID() << '\t' <<  elems[ee]->getLevel() << endl;

           elems[ee]->resetMatrixAndVector();
           elems[ee]->setInitialProfile();
           //elems[ee]->AssembleMatrixAndVector(1, solver->mtx, &(rhsVec(0)));
       }
    }

    solverEigen->currentStatus = ASSEMBLY_OK;

    factoriseSolveAndUpdate();

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



int HBSplineCutFEM::calcStiffnessAndResidual(int solver_type, bool zeroMtx, bool zeroRes)
{
    PetscPrintf(MPI_COMM_WORLD, "     HBSplineCutFEM: generating coefficient Matrices ...\n\n");

    assert( SOLVER_TYPE == SOLVER_TYPE_PETSC );

    //char fct[] = "HBSplineCutFEM::calcStiffnessAndResidual";
    //computerTime.go(fct);

    int  ii, ee, nr, nc, dd, start=0, pp, dof, dir, domTemp=0, bb, kk;
    double  PENALTY;

    IterNum   = (iterCount == 1);

    if(firstIter)
    {
      //solver->ResetPrecondFlag();
      rNorm = -1.0;
    }

    ////////////////////////////////////////////////////////
    // stiffness and residual for the background fluid grid
    ////////////////////////////////////////////////////////

    //tstart = time(0);
    //tstart = MPI_Wtime();
    auto tstart = Clock::now();

    solverPetsc->zeroMtx();

    //cout << " wwwwwwwwwww " << endl;

    //#pragma omp parallel for
    for(ee=0; ee<fluidElementIds.size(); ee++)
    {
      node *nd = elems[fluidElementIds[ee]];

      //cout << " nd->getID() " <<  nd->getID() << '\t' <<  nd->getLevel() << '\t' << nd->getDomainNumber() << endl;

      if( nd->getSubdomainId() == this_mpi_proc )
      {
        domTemp = nd->getDomainNumber() ;

        if( domTemp <= 0 )
        {
          nr = nd->forAssyVec.size();

          //cout << nr << '\t' << nc << endl;
          MatrixXd  Klocal;
          VectorXd  Flocal;

          Klocal = MatrixXd::Zero(nr, nr);
          Flocal = VectorXd::Zero(nr);

          //cout << " AAAAAAAAAAAAAAAAA " << endl;
          //nd->calcStiffnessAndResidualGFEM(Klocal, Flocal, domTemp);
          nd->calcStiffnessAndResidualCutFEMFluid(Klocal, Flocal, domTemp);

          //cout << " BBBBBBBBBBBBBBBBB " << endl;
          //nd->applyDirichletBCsGFEM(Klocal, Flocal, domTemp);
          nd->applyDirichletBCsCutFEMFluid(Klocal, Flocal, domTemp);
          //cout << " DDDDDDDDDDDDDDDDD " << endl;

          //cout << " DDDDDDDDDDDDDDDDD " << endl;
          nd->applyNeumannBCsCutFEMFluid(Klocal, Flocal, domTemp);
          //cout << " BBBBBBBBBBBBBBBBB " << endl;

          //printMatrix(Klocal);
          //printf("\n\n\n");
          //printVector(Flocal);

          //#pragma omp critical
          solverPetsc->AssembleMatrixAndVectorCutFEM(start, start, nd->forAssyVec, grid_to_proc_DOF, Klocal, Flocal);
          //cout << " CCCCCCCCCCCCCCCC " << endl;
        }
      }
    }

    //cout << " MMMMMMMMMMMMMMMM " << endl;
    //printVector(pointBCs[0]);

    myDataIntegrateCutFEM  myData;

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

        geometryToParametric(geom, param);
        
        myData.param = param;
        myData.geom  = geom;

        node *nd;
        nd = elems[ee];

        int nr = nd->forAssyVec.size();

        //MatrixXd  Klocal;
        //VectorXd  Flocal;

        //Klocal = MatrixXd::Zero(nr, nc);
        //Flocal = VectorXd::Zero(nr);

        myData.K1 = MatrixXd::Zero(nr, nr);
        myData.F1 = VectorXd::Zero(nr);

        //dof, param, spec_val, PENALTY, Klocal, Flocal
        nd->applyBoundaryConditionsAtApoint(myData);
        solverPetsc->AssembleMatrixAndVector(0, 0, nd->forAssyVec, nd->forAssyVec, myData.K1, myData.F1);
      }
    }


    //cout << " rhsVec " << endl;
    //       printVector(&(solver->rhsVec(0)), totalDOF);

    //printf("\n rhsVec norm = %12.6E \n", solverEigen->rhsVec.norm());

    if(DIM == 2)
      applyInterfaceTerms2D();

    if(DIM == 3)
      applyInterfaceTerms3D();

    //printf("\n rhsVec norm = %12.6E \n", solverEigen->rhsVec.norm());

    if(!STAGGERED)
    {
      kk = fluidDOF;
      for(bb=0; bb<ImmersedBodyObjects.size(); bb++)
      {
        cout << " ppppppppppp " << kk << endl;
        ImmersedBodyObjects[bb]->AssembleGlobalMatrixAndVectorCutFEM(kk, kk, solverPetsc);
        cout << " ppppppppppp " << endl;
        kk += ImmersedBodyObjects[bb]->GetTotalDOF();
      }
    }

    //cout << " rhsVec " << endl;
    //printVector(&(solver->rhsVec(0)), totalDOF);
    //printf("\n rhsVec norm = %12.6E \n", solver->rhsVec.norm());

    firstIter = false;
    rNormPrev = rNorm;

    VecAssemblyBegin(solverPetsc->rhsVec);
    VecAssemblyEnd(solverPetsc->rhsVec);

    VecNorm(solverPetsc->rhsVec, NORM_2, &rNorm);
    solverPetsc->currentStatus = ASSEMBLY_OK;

    //IterNum   = (iterCount == 1);

    SolnData.firstIter = firstIter;

    if( std::isnan(rNorm) )
    {
      //VecView(solverPetsc->rhsVec, PETSC_VIEWER_STDOUT_WORLD);

      //printVector(grid_to_proc_DOF);
      //printVector(SolnData.var1);

      cerr << "  NAN found "  << endl;
      exit(0);
    }

    PetscPrintf(MPI_COMM_WORLD, "  %5d \t %11.4e\n", iterCount, rNorm);

    //ctimCalcStiffRes += computerTime.stop(fct);
    //computerTime.stopAndPrint(fct);

    //PetscPrintf(MPI_COMM_WORLD, "     HBSplineCutFEM: generating coefficient Matrices ...DONE  \n\n");

    auto tend = Clock::now();
    //tend = time(0);
    //tend = MPI_Wtime();
    PetscPrintf(MPI_COMM_WORLD, "HBSplineCutFEM::calcStiffnessAndResidual() took %d millisecond(s) \n\n ", std::chrono::duration_cast<std::chrono::milliseconds>(tend - tstart).count());

    iterCount++;

    return 0;
}






int HBSplineCutFEM::factoriseSolveAndUpdate()
{
  //cout << "     HBSplineCutFEM: solving the matrix system ...  \n\n";

  int  bb, ii, jj, dd, kk=0;

  //solver->computeConditionNumber();

  //cout << " rhsVec " << endl;        printVector(&(solver->rhsVec(0)), totalDOF);

  //auto tstart = Clock::now();
  Clock::time_point tstart = Clock::now();


  //VecView(solverPetsc->rhsVec, PETSC_VIEWER_STDOUT_WORLD);

  //tstart = time(0);

  solverPetsc->factoriseAndSolve();

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

  kk=0;
  for(ii=0; ii<nNode; ii++)
  {
    //cout << ii << '\t' << node_map_old_to_new[ii] << '\t' << node_map_new_to_old[ii] << endl;
    //jj = node_map_old_to_new[ii]*ndof;
    jj = ii*ndof;
    for(dd=0; dd<ndof; dd++)
      soln[kk++] = arrayTemp[jj+dd];
  }

  for(ii=0; ii<fluidDOF; ii++)
    SolnData.var1[proc_to_grid_DOF[ii]] += soln[ii];

    if(!STAGGERED)
    {
      cout << " need to update this for monolithic scheme " << endl;
      kk = fluidDOF;
      for(bb=0; bb<ImmersedBodyObjects.size(); bb++)
      {
        //if(ImmersedBodyObjects[bb]->getNdof() > 0)
        //{
          cout << " AAAAAAAAAAA " << bb << endl;
          ImmersedBodyObjects[bb]->updateDisplacement(&arrayTemp[kk]);
          //cout << " AAAAAAAAAAA " << bb << endl;
          kk += ImmersedBodyObjects[bb]->GetTotalDOF();
        //}
      }
    }

  VecRestoreArray(vec_SEQ, &arrayTemp);

  //printVector(SolnData.var1);
  //printf("\n\n\n");

  //computerTime.stopAndPrint(fct);
  //ctimFactSolvUpdt += computerTime.stop(fct);

  //cout << " result " << endl;        printVector(&(soln(0)), totalDOF);

  //cout << "     HBSplineCutFEM: solving the matrix system ...DONE  \n\n";

  //auto tend = Clock::now();
  Clock::time_point tend = Clock::now();
  //tend = MPI_Wtime();
  //cout << " HBSplineCutFEM::factoriseSolveAndUpdate() took " << std::chrono::duration_cast<std::chrono::milliseconds>(tend - tstart).count() << " millisecond(s) " << endl;
  PetscPrintf(MPI_COMM_WORLD, "HBSplineCutFEM::factoriseSolveAndUpdate() took %d millisecond(s) \n ", std::chrono::duration_cast<std::chrono::milliseconds>(tend - tstart).count());

  //computerTime.stopAndPrint(fct);

  return 0;
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
            domTemp = nd->getDomainNumber() ;
            //if(domTemp == 1)
            //{
              nd->calcError(index, domTemp);
            
              //cout << ee << '\t' << elems[ee]->getError() << endl;
              //totalError += ( elems[ee]->getError() * elems[ee]->getError() );
              totalError +=  nd->getError();
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
            
          totalError += nd->getError();

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





