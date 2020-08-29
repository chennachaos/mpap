
#include "HBSplineCutFEM.h"
#include "SolverEigen.h"
#include "ComputerTime.h"
#include "MpapTime.h"
#include "myDataIntegrateCutFEM.h"
#include "ImmersedIntegrationElement.h"
#include "TimeFunction.h"


extern ComputerTime       computerTime;
extern MpapTime mpapTime;
extern List<TimeFunction> timeFunction;



void HBSplineCutFEM::setInitialConditions()
{
    PetscPrintf(MPI_COMM_WORLD, "     HBSplineCutFEM: setting initial conditions ...\n\n");

    SolnData.var1.setZero();

    VectorXd  specVal(3);
    specVal.setZero();
    specVal(0) = DirichletBCs[0][2] * timeFunction[0].prop;
    specVal(2) = 0.1;

    //timeUpdate();

    solverPetsc->zeroMtx();

    int start = 0, nr, domTemp;
    for(int ee=0; ee<fluidElementIds.size(); ee++)
    {
      node *nd = elems[fluidElementIds[ee]];

      if( nd->getSubdomainId() == this_mpi_proc )
      {
        domTemp = nd->getDomainNumber() ;

        if( domTemp <= 0 )
        {
          nr = nd->forAssyVec.size();

          MatrixXd  Klocal = MatrixXd::Zero(nr, nr);
          VectorXd  Flocal = VectorXd::Zero(nr);

          nd->setInitialProfile(Klocal, Flocal, specVal, domTemp);

          solverPetsc->assembleMatrixAndVectorCutFEM(start, start, nd->forAssyVec, grid_to_proc_DOF, Klocal, Flocal);
        }
      }
    }

    if(DIM == 2)
      applyInterfaceTerms2D();

    solverPetsc->currentStatus = ASSEMBLY_OK;

    PetscPrintf(MPI_COMM_WORLD, "     HBSplineCutFEM: solving initial conditions ...\n\n");

    factoriseSolveAndUpdate();

    //for(int ii=0; ii<gridBF1; ii++)
      //SolnData.var1[ii*ndof] = DirichletBCs[0][2];

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
    assert( SOLVER_TYPE == SOLVER_TYPE_PETSC );

    //PetscPrintf(MPI_COMM_WORLD, "     HBSplineCutFEM: generating coefficient Matrices ...\n\n");

    /*
    PetscLogDouble mem1, mem2, mem3, mem4;
    //the current resident set size (memory used) for the program.
    ierr = PetscMemoryGetCurrentUsage(&mem1);               //CHKERRQ(ierr);
    //the maximum resident set size (memory used) for the program.
    ierr = PetscMemoryGetMaximumUsage(&mem2);               //CHKERRQ(ierr);
    //the current amount of memory used that was PetscMalloc()ed
    ierr = PetscMallocGetCurrentUsage(&mem3);               //CHKERRQ(ierr);
    //the maximum amount of memory used that was PetscMalloc()ed at any time during this run.
    ierr = PetscMallocGetMaximumUsage(&mem4);               //CHKERRQ(ierr);

    //PetscPrintf(MPI_COMM_WORLD, " Petsc memory allocation details ... %12.8f \t %12.8f \t%12.8f \t%12.8f \n\n", mem1, mem2, mem3, mem4);
    */

    int  ii=0, ee=0, nr=0, nc=0, dd=0, start=0, pp=0;
    int  dof=0, dir=0, domTemp=0, bb=0, kk=0;
    double  PENALTY=0.0;

    IterNum   = (iterCount == 1);

    if(firstIter)  rNorm = -1.0;

    ////////////////////////////////////////////////////////
    // stiffness and residual for the background fluid grid
    ////////////////////////////////////////////////////////

    double tstart = MPI_Wtime();

    // create one matrix and vector of large enough size 
    // to compute element stiffness and residual
    // This way, we avoid storing, or creating a new one,
    // local matrix/vector for every element
    //MatrixXd  Klocal(50,50);
    //VectorXd  Flocal(50);

    solverPetsc->zeroMtx();

    for(ee=0; ee<fluidElementIds.size(); ee++)
    {
      node *nd = elems[fluidElementIds[ee]];

      if( nd->getSubdomainId() == this_mpi_proc )
      {
        domTemp = nd->getDomainNumber() ;

        if( domTemp <= 0 )
        {
          nr = nd->forAssyVec.size();

          MatrixXd  Klocal = MatrixXd::Zero(nr, nr);
          VectorXd  Flocal = VectorXd::Zero(nr);
          //Klocal = MatrixXd::Zero(nr, nr);
          //Flocal = VectorXd::Zero(nr);
          //Klocal.setZero();
          //Flocal.setZero();

          nd->calcStiffnessAndResidualCutFEMFluid(Klocal, Flocal, domTemp);

          nd->applyDirichletBCsCutFEMFluid(Klocal, Flocal, domTemp);

          nd->applyNeumannBCsCutFEMFluid(Klocal, Flocal, domTemp);

          nd->applyDerivativeBCsCutFEMFluid(Klocal, Flocal, domTemp);

          solverPetsc->assembleMatrixAndVectorCutFEM(start, start, nd->forAssyVec, grid_to_proc_DOF, Klocal, Flocal);
        }
      }
    }

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

        nr = nd->forAssyVec.size();

        myData.K1 = MatrixXd::Zero(nr, nr);
        myData.F1 = VectorXd::Zero(nr);

        //dof, param, spec_val, PENALTY, Klocal, Flocal
        nd->applyBoundaryConditionsAtApoint(myData);
        solverPetsc->assembleMatrixAndVector(0, 0, nd->forAssyVec, nd->forAssyVec, myData.K1, myData.F1);
      }
    }

    //cout << " rhsVec " << endl;
    //printVector(&(solver->rhsVec(0)), totalDOF);
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
        //cout << " ppppppppppp " << kk << endl;
        ImmersedBodyObjects[bb]->assembleGlobalMatrixAndVectorCutFEM(kk, kk, solverPetsc);
        //cout << " ppppppppppp " << endl;
        kk += ImmersedBodyObjects[bb]->getTotalDOF();
      }
    }

    //printVector(&(solver->rhsVec(0)), totalDOF);
    //printf("\n rhsVec norm = %12.6E \n", solverPetsc->rhsVec.norm());

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
      VecView(solverPetsc->rhsVec, PETSC_VIEWER_STDOUT_WORLD);

      cerr << "  NAN found "  << endl;
      exit(0);
    }

    PetscPrintf(MPI_COMM_WORLD, "HBSplineCutFEM-- %d \t %11.4e\n", iterCount, rNorm);

    if(IB_MOVED)
    {
      for(ee=0; ee<cutCellIds.size(); ee++)
      {
        elems[cutCellIds[ee]]->Quadrature.reset();
      }
    }

    double tend = MPI_Wtime();
    PetscPrintf(MPI_COMM_WORLD, " HBSplineCutFEM::calcStiffnessAndResidual() took %f  milliseconds \n", (tend-tstart)*1000);

    iterCount++;

    return 0;
}






int HBSplineCutFEM::factoriseSolveAndUpdate()
{
    //cout << "     HBSplineCutFEM: solving the matrix system ...  \n\n";

    int  bb, ii, jj, dd, kk=0;

    //cout << " rhsVec " << endl;        printVector(&(solver->rhsVec(0)), totalDOF);

    //VecView(solverPetsc->rhsVec, PETSC_VIEWER_STDOUT_WORLD);

    double tstart = MPI_Wtime();

    if( solverPetsc->factoriseAndSolve() )
    {
        int  resln[] = {1, 1, 1};

        postProcessFlow(1, 1, 1, 1, 0.0, 1.0, resln);

        return 1;
    }

    double tend = MPI_Wtime();
    PetscPrintf(MPI_COMM_WORLD, " PETSC solver took %f  milliseconds \n", (tend-tstart)*1000);

    //VecView(solver2->soln, PETSC_VIEWER_STDOUT_WORLD);

    /////////////////////////////////////////////////////////////////////////////
    // get the solution vector onto all the processors
    /////////////////////////////////////////////////////////////////////////////

    Vec            vec_SEQ;
    VecScatter     ctx;
    PetscScalar *arrayTemp;

    if(n_mpi_procs > 1)
    {
      VecScatterCreateToAll(solverPetsc->soln, &ctx, &vec_SEQ);
      VecScatterBegin(ctx, solverPetsc->soln, vec_SEQ, INSERT_VALUES, SCATTER_FORWARD);
      VecScatterEnd(ctx, solverPetsc->soln, vec_SEQ, INSERT_VALUES, SCATTER_FORWARD);
      VecGetArray(vec_SEQ, &arrayTemp);
    }
    else
    {
      VecGetArray(solverPetsc->soln, &arrayTemp);
    }

    // update solution vector

    kk=0;
    for(ii=0; ii<nNode; ii++)
    {
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
          kk += ImmersedBodyObjects[bb]->getTotalDOF();
        //}
      }
    }

    if(n_mpi_procs > 1)
    {
      VecRestoreArray(vec_SEQ, &arrayTemp);
      VecScatterDestroy(&ctx);
      VecDestroy(&vec_SEQ);
    }
    else
    {
      VecRestoreArray(solverPetsc->soln, &arrayTemp);
    }

    //double tend = MPI_Wtime();
    //PetscPrintf(MPI_COMM_WORLD, " HBSplineCutFEM::factoriseSolveAndUpdate() took %f  milliseconds \n", (tend-tstart)*1000);

    return 0;
}




void HBSplineCutFEM::computeElementErrors(int index)
{
    if(index==5)
    {
      solverEigen->computeConditionNumber();
      return;
    }

    int  ii=0, ee=0, count=0, dd=0, domTemp=0;
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

              //cout << ee << '\t' << elems[ee]->getError() << endl;
              totalError +=  nd->calcError(index, domTemp);
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

          totalError += nd->calcError(index);

          count++;
       }
       totalError /= count;
    }

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





