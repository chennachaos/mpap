

#include "HBSplineFEM.h"
#include "ComputerTime.h"
#include "MpapTime.h"
#include "myDataIntegrateCutFEM.h"
#include "ContactElementPointToPoint2D.h"
#include "QuadratureUtil.h"
#include "headersEigen.h"


#include <omp.h>


extern ComputerTime       computerTime;
extern MpapTime mpapTime;


void HBSplineFEM::setInitialConditions()
{
    double* tmp;

    solverEigen->zeroMtx();

    for(int ii=0;ii<activeElements.size();ii++)
    {
       //cout << ii << '\t' << omp_get_thread_num() << '\t' << omp_get_num_threads() << '\t' << omp_get_max_threads() << endl;
       node *nd = elems[activeElements[ii]];

       nd->setInitialProfile();

       //cout << " AAAAAAAAAAAAAAAAA " << endl;

       //cout << " BBBBBBBBBBBBBBBBB " << endl;
       //printMatrix(Klocal);
       //printf("\n\n");
       //printVector(Flocal);

      //nd->assembleElementMatrix2(0, solverEigen->mtx);
      //nd->assembleElementVector(0, 0, &(solverEigen->rhsVec(0)));
    }

    solverEigen->currentStatus = ASSEMBLY_OK;

    factoriseSolveAndUpdate();

    //
    int resln1[3]; resln1[0] = resln1[1] = 5;

    //postProcess2D(1, 1, 10, 1, 0.0, 1.0, resln1);
    postProcess1D(1, 1, 10, 1, 0.0, 1.0, resln1);
    //

    //printVector(solver->soln);
    
    //solnInit = soln;
    //SolnData.var1 = solnInit;
    //SolnData.var1Prev  = SolnData.var1;
    //SolnData.var1Prev2 = SolnData.var1;
    //SolnData.var1Prev3 = SolnData.var1;
    //SolnData.var1Prev4 = SolnData.var1;

    //printVector(SolnData.var1);

   return;
}


/*
int HBSplineFEM::calcStiffnessAndResidual(int solver_type, bool zeroMtx, bool zeroRes)
{
    ////////////////////////////////////////////////////////
    //
    // stiffness and residual for the LSFEM
    //
    ////////////////////////////////////////////////////////

    //cout << "     HBSplineFEM: generating coefficient Matrices ...\n\n";

    char fct[] = "HBSplineFEM::calcStiffnessAndResidual";
    computerTime.go(fct);

    int  ee, ii, jj, kk, ll, nr, nc, aa, bb, dd, numPoints, size2, r, c;
    Point  val(3);
    time_t tstart, tend;
    
    //MatrixXd  Klocal;
    //VectorXd  Flocal;

    if(firstIter)
       rNorm = -1.0;

    solver->zeroMtx();
    
    ////////////////////////////////////////////////////////
    //
    // stiffness and residual for the background fluid grid
    //
    ////////////////////////////////////////////////////////

    tstart = time(0);
    //node *nd;

    //printf("Number of threads: %i \n",omp_get_num_procs());
    //omp_set_num_threads(1);
    //#pragma omp parallel
    //{
      //printf("Number of threads: %i \n",omp_get_num_procs());
    //}


  if(LSFEM_FLAG)
  {
    //#pragma omp parallel for
    for(ii=0;ii<activeElements.size();ii++)
    {
       //cout << ii << '\t' << omp_get_thread_num() << '\t' << omp_get_num_threads() << '\t' << omp_get_max_threads() << endl;
       node *nd;
       nd = elems[activeElements[ii]];
       //cout << " nd->getID() " <<  nd->getID() << '\t' <<  nd->getLevel() << endl;

       //nd->resetMatrixAndVector();
       //FluidSolnData.resetMatrixAndVector(nd->getNsize2());
       int nr, nc;

       //printVector(nd->forAssyVec);
       //printVector(nd->forAssyVec2);
       
       nr = nd->forAssyVec.size();
       nc = nr;

       //cout << nr << '\t' << nc << endl;
       MatrixXd  Klocal;
       VectorXd  Flocal;

       Klocal = MatrixXd::Zero(nr, nc);
       Flocal = VectorXd::Zero(nr);
       //Klocal.resize(nr, nc);
       //Flocal.resize(nr);
       //Klocal.setZero();
       //Flocal.setZero();

       //cout << " AAAAAAAAAAAAAAAAA " << endl;
       nd->calcStiffnessAndResidualLSFEM(Klocal, Flocal);
       //cout << " AAAAAAAAAAAAAAAAA " << endl;
       //printMatrix(Klocal);
       //printf("\n\n");
       //printVector(Flocal);
       nd->applyDirichletBCsLSFEM(Klocal, Flocal);
       nd->applyNeumannBCsLSFEM(Klocal, Flocal);
       //cout << " AAAAAAAAAAAAAAAAA " << endl;
       //printMatrix(Klocal);
       //printf("\n\n");
       //printVector(Flocal);
       
       //nd->assembleMatrixAndVector(1, ((SolverEigen*)solver)->mtx, &(rhsVec(0)));
       //nd->assembleElementMatrix(1, ((SolverEigen*)solver)->mtx);
       //nd->assembleElementVector(firstIter, 1, &(rhsVec(0)));
       //cout << " AAAAAAAAAAAAAAAAA " << endl;
       //#pragma omp critical
        //solver->assembleMatrixAndVector(nd->forAssyVec, nd->forAssyVec, Klocal, Flocal);
        solver->assembleMatrixAndVector(velDOF, 0, nd->forAssyVec, Klocal, Flocal);
       //cout << " AAAAAAAAAAAAAAAAA " << endl;
    }
  }

    //cout << " rhsVec " << endl;        printVector(&(solver->rhsVec(0)), totalDOF);

    //applyBoundaryConditions();
    //addExternalForces();

    tend = time(0); 
    cout << "HBSplineFEM::calcStiffnessAndResidual()  took "<< difftime(tend, tstart) <<" second(s)."<< endl;

    printf("\n rhsVec norm = %12.6E \n", solver->rhsVec.norm());

    //printData(1, 1);
    //printData(3, 1);

    ////////////////////////////////////////////////////////
    //
    // stiffness and residual for the immersed boundary points
    // Lagrange multipliers
    //  or
    // Penalty method
    //
    ////////////////////////////////////////////////////////

    //cout << " AAAAAAAAAAAAAAAAA " << endl;
    //cout << " AAAAAAAAAAAAAAAAA " << endl;

    //ImmersedBoundaryBodyForce();

    MatrixXd  Klocal;
    VectorXd  Flocal;
    node *nd;
    Bpoint *bp;

    if(!LSFEM_FLAG)
    {
      for(bb=0;bb<ImmersedBodyObjects.size();bb++)
      {
        if(ImmersedBodyObjects[bb]->isBoundaryConditionTypeLagrange())
        {
          for(aa=0;aa<ImmersedBodyObjects[bb]->getNumberOfNodes();aa++)
          {
            //cout << bb << '\t' << aa << '\t' << ii << endl;

            //nr = nd->getNsize2();
            //nc = nr;

            //Klocal = MatrixXd::Zero(nr, nc);
            //Flocal = VectorXd::Zero(nr);

            ImmersedBodyObjects[bb]->IBpoints[aa]->resetMatrixAndVector();
            //cout << " PPPPPPPPPPPPPPPPPPPP " << endl;
            ImmersedBodyObjects[bb]->IBpoints[aa]->calcStiffnessAndResidual(1,0,0.0,0.0);
            //cout << " PPPPPPPPPPPPPPPPPPPP " << endl;
            ImmersedBodyObjects[bb]->IBpoints[aa]->assembleMatrixAndVector(fluidDOF, solver->mtx, &(solver->rhsVec(0)));
            //solver->assembleMatrixAndVector(nd->forAssyVec, nd->forAssyVec, Klocal, Flocal);
          }
        }
        else
        {
          for(aa=0;aa<ImmersedBodyObjects[bb]->IBpoints.size();aa++)
          {
            bp = ImmersedBodyObjects[bb]->IBpoints[aa];
            //cout << " uuuuuuuuuuu " << endl;
            ee = bp->GetElementNum();
            nd = elems[bp->GetElementNum()];

            for(dd=0;dd<DIM;dd++)
            {
              param[dd] = bp->GetParam(dd);
              val[dd]   = bp->GetSpecVal(dd);
            }

            //cout << aa << '\t' << ee << '\t' << param[0] << '\t' << param[1] << '\t' << val[0] << '\t' << val[1] << endl;

            nr = nd->getNsize2();
            nc = nr;

            Klocal = MatrixXd::Zero(nr, nc);
            Flocal = VectorXd::Zero(nr);

            //cout << " uuuuuuuuuuu " << endl;
            for(jj=0;jj<DIM;jj++)
              nd->applyBoundaryConditionsAtApoint(jj, param, val[jj], Klocal, Flocal);

            solver->assembleMatrixAndVector(nd->forAssyVec, nd->forAssyVec, Klocal, Flocal);
            //nd->assembleMatrixAndVector(1, ((SolverEigen*)solver)->mtx, &(rhsVec(0)));
            //cout << " uuuuuuuuuuu " << endl;
          }
        }
      }
    }

    //printData(1, 1);
    //printData(3, 1);

    printf("\n rhsVec norm = %12.6E \n", solver->rhsVec.norm());

   //cout << " rhsVec " << endl;
   //for(int ii=fluidDOF-10;ii<totalDOF;ii++)
     //printf("%5d \t %12.8f \n", ii, rhsVec(ii));

    //printf("\n rhsVec norm = %12.6E \n", solver->rhsVec.norm());

    //cout << " rhsVec " << endl;
    //for(int ii=fluidDOF-10;ii<totalDOF;ii++)
      //printf("%5d \t %12.8f \n", ii, rhsVec(ii));

    firstIter = false;
    rNormPrev = rNorm;
    rNorm     = solver->rhsVec.norm();
    iterCount++;
    
    COUT << domain.name(this); printf("  %11.4e\n",rNorm);

    if(IterNum)
      solver->currentStatus = ASSEMBLY_OK; //FACTORISE_OK

    ctimCalcStiffRes += computerTime.stop(fct);
    //computerTime.stopAndPrint(fct);
   
    //cout << "     HBSplineFEM: generating coefficient Matrices ...DONE  \n\n";
   
    return 0;
}
*/



int HBSplineFEM::calcStiffnessAndResidual(int solver_type, bool zeroMtx, bool zeroRes)
{
    //cout << "     HBSplineFEM: generating coefficient Matrices ...\n\n";

    char fct[] = "HBSplineFEM::calcStiffnessAndResidual";
    computerTime.go(fct);

    int bb, ee, ii, jj, kk, ll, nr, nc, aa, dd, numPoints, size2, r, c, gp;

    time_t tstart, tend;
    
    IterNum   = (iterCount == 1);

    //cout << " IterNum " << IterNum << endl;

    if(firstIter)
    {
      rNorm = -1.0;
      //solver->mtx *= 0.0;
    }

    solverEigen->zeroMtx();

    ////////////////////////////////////////////////////////
    //
    // stiffness and residual for the background fluid grid
    //
    ////////////////////////////////////////////////////////

    //#pragma omp parallel
      //printf("Number of threads: %i \n",omp_get_num_procs());

      MatrixXd  matM, matK;

    tstart = time(0);
    //#pragma omp parallel for
    for(ii=0;ii<activeElements.size();ii++)
    {
      //cout << ii << '\t' << omp_get_thread_num() << '\t' << omp_get_num_threads() << '\t' << omp_get_max_threads() << endl;
      node *nd = elems[activeElements[ii]];
      //cout << " nd->getID() " <<  nd->getID() << '\t' <<  nd->getLevel() << endl;

      //printVector(nd->forAssyVec);
      //printVector(nd->forAssyVec2);

      int nr = nd->forAssyVec.size();

      //cout << nr << '\t' << nc << endl;
      MatrixXd  Klocal;
      VectorXd  Flocal;

      Klocal = MatrixXd::Zero(nr, nr);
      Flocal = VectorXd::Zero(nr);

      matM = MatrixXd::Zero(nr, nr);
      matK = MatrixXd::Zero(nr, nr);

      //cout << " AAAAAAAAAAAAAAAAA " << endl;
      nd->calcStiffnessAndResidualGFEM(Klocal, Flocal);
      matK = Klocal;
      //cout << " AAAAAAAAAAAAAAAAA " << endl;
      //printMatrix(Klocal);
      //printf("\n\n");
      //printVector(Flocal);
      Klocal.setZero();
      nd->applyDirichletBCsGFEM(Klocal, Flocal);
      matM = Klocal;
      //cout << " BBBBBBBBBBBBBBBBB " << endl;
      //nd->applyNeumannBCsGFEM(Klocal, Flocal);
      //cout << " BBBBBBBBBBBBBBBBB " << endl;
      //printMatrix(Klocal);
      //printf("\n\n");
      //printVector(Flocal);

       //cout << " AAAAAAAAAAAAAAAAA " << endl;
       //nd->calcStiffnessAndResidualLSFEM(1, Klocal, Flocal);
       //cout << " AAAAAAAAAAAAAAAAA " << endl;
       //printMatrix(Klocal);
       //printf("\n\n");
       //printVector(Flocal);
       //nd->applyDirichletBCsLSFEM(1, Klocal, Flocal);
       //nd->applyNeumannBCsLSFEM(1, Klocal, Flocal);
       //cout << " AAAAAAAAAAAAAAAAA " << endl;

       //solver->assembleMatrixAndVector(velDOF, 0, nd->forAssyVec, nd->forAssyVec2, Klocal, Flocal);

       //if(IterNum)
       //cout << " AAAAAAAAAAAAAAAAA " << endl;
      //#pragma omp critical
        solverEigen->assembleMatrixAndVector(0, 0, nd->forAssyVec, Klocal, Flocal);
      //else
        //solverEigen->assembleVector(velDOF, 0, nd->forAssyVec, Flocal);
      //cout << " AAAAAAAAAAAAAAAAA " << endl;
    }

    //solverEigen->computeConditionNumber();
    //cout << " AAAAAAAAAAAAAAAAA " << endl;
    // for computing CI
    //GeneralizedSelfAdjointEigenSolver<MatrixXd> es(matK, matM);
    // for computing penalty parameter
    GeneralizedSelfAdjointEigenSolver<MatrixXd> es(matM, matK);
    //cout << " AAAAAAAAAAAAAAAAA " << endl;
    VectorXd  eig_vals = es.eigenvalues();
    cout << "The eigenvalues of the pencil (A,B) are:" << endl << eig_vals << endl;

    printf("\n %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");
    printf("\n\n\n Maximum eigenvalue = %12.6f \n\n\n", eig_vals.maxCoeff() );
    //printf("\n\n\n CI value = %14.10f \n\n\n", pow(eig_vals.maxCoeff(),1.0/3.0) );
    printf("\n %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");

    //
    if(pointBCs.size() > 0)
    {
      double  spec_val, PENALTY;
      int ii, pp, ee, dof, dir;
    
      for(pp=0;pp<pointBCs.size();pp++)
      {
        geom[0] = pointBCs[pp][0];
        geom[1] = pointBCs[pp][1];
        geom[2] = pointBCs[pp][2];

        dof = int (pointBCs[pp][3] - 1);

        spec_val = pointBCs[pp][4];
        PENALTY  = pointBCs[pp][5];
        //printVector(pointBCs[0]);

        //cout << findCellNumber(geom) << endl;

        ee = findCellNumber(geom);

        geometryToParametric(geom, param);

        node *nd;
        nd = elems[ee];

        int nr = nd->forAssyVec.size();// + nd->forAssyVec2.size();
        int nc = nr;

        //cout << nr << '\t' << nc << '\t' << spec_val << '\t' << PENALTY << endl;
        MatrixXd  Klocal;
        VectorXd  Flocal;

        Klocal = MatrixXd::Zero(nr, nc);
        Flocal = VectorXd::Zero(nr);

        cout << " this needs to completed " << endl;
        //nd->applyBoundaryConditionsAtApoint(dof, param, spec_val, PENALTY, Klocal, Flocal);
        solverEigen->assembleMatrixAndVector(velDOF, 0, nd->forAssyVec, Klocal, Flocal);
      }
    }
    //
    //cout << " rhsVec " << endl;        printVector(&(solver->rhsVec(0)), totalDOF);

    //applyBoundaryConditions();
    //addExternalForces();

    tend = time(0); 
    //cout << "HBSplineFEM::calcStiffnessAndResidual()  took "<< difftime(tend, tstart) <<" second(s)."<< endl;

    printf("\n rhsVec norm = %12.6E \n", solverEigen->rhsVec.norm());

    ////////////////////////////////////////////////////////
    //
    // stiffness and residual for the immersed boundary points
    // Lagrange multipliers
    //  or
    // Penalty method
    //
    ////////////////////////////////////////////////////////

    //cout << " AAAAAAAAAAAAAAAAA " << velDOF << '\t' << presDOF << '\t' << fluidDOF << endl;

    if(DIM == 2)
      applyInterfaceTerms2D();

    if(DIM == 3)
      applyInterfaceTerms3D();

    //printData(1, 1);
    //printData(3, 1);
    //cout << " rhsVec " << endl;        printVector(&(solver->rhsVec(0)), totalDOF);
    printf("\n rhsVec norm = %12.6E \n", solverEigen->rhsVec.norm());

    //cout << " rhsVec " << endl;
   //for(int ii=fluidDOF-10;ii<totalDOF;ii++)
     //printf("%5d \t %12.8f \n", ii, rhsVec(ii));

    ////////////////////////////////////////////////////////
    //
    // stiffness and residual for the immersed solid body
    // only for the MONOLITHIC_SCHEME
    //
    ////////////////////////////////////////////////////////

    //cout << " uuuuuuuuuuu " << endl;
    if(!STAGGERED)
    {
      //cout << " Terms related to monolithic scheme " << endl;

      for(bb=0;bb<ImmersedBodyObjects.size();bb++)
      {
          kk = fluidDOF + IBDOF;
          //cout << " ppppppppppp " << kk << endl;
          ImmersedBodyObjects[bb]->assembleGlobalMatrixAndVector(fluidDOF, kk, solverEigen->mtx, &(solverEigen->rhsVec(0)));
          //cout << " ppppppppppp " << endl;
      }

      kk = fluidDOF + IBDOF;
      for(bb=0;bb<contactElementObjects.size();bb++)
      {
        contactElementObjects[bb]->calcStiffnessAndResidual();
        contactElementObjects[bb]->assembleMatrixAndVector(kk, solverEigen->mtx, &(solverEigen->rhsVec(0)));
      }
    }

    //printf("\n rhsVec norm = %12.6E \n", solverEigen->rhsVec.norm());

    //cout << " rhsVec " << endl;
    //for(int ii=fluidDOF-10;ii<totalDOF;ii++)
      //printf("%5d \t %12.8f \n", ii, rhsVec(ii));

    firstIter = false;
    rNormPrev = rNorm;
    rNorm     = solverEigen->rhsVec.norm();
    //IterNum   = (iterCount == 1);
    iterCount++;
    
    SolnData.firstIter = firstIter;

    if(IBDOF > 0)
    {
      if(!STAGGERED)
      {
        if(ImmersedBodyObjects.size() > 0 )
          for(int bb=0;bb<ImmersedBodyObjects.size();bb++)
            ImmersedBodyObjects[bb]->firstIter = firstIter;
      }
    }

    printf(" %5d \t %11.4e\n", iterCount, rNorm);
    //COUT << domain.name(this); printf("  %11.4e\n",rNorm);

    //if(IterNum)
      solverEigen->currentStatus = ASSEMBLY_OK;
    //else
      //solverEigen->currentStatus = FACTORISE_OK;

    ctimCalcStiffRes += computerTime.stop(fct);
    //computerTime.stopAndPrint(fct);

    //cout << "     HBSplineFEM: generating coefficient Matrices ...DONE  \n\n";

    return 0;
}





void  HBSplineFEM::applyBoundaryConditions()
{
    //cout << "     HBSplineFEM: applyBoundaryConditions ... STARTED  \n\n";
    //firstIter = true;

    //cout << "     HBSplineFEM: applyBoundaryConditions ...DONE  \n\n";
    return;
}



void HBSplineFEM::addExternalForces()
{
//  cout << "HBSplineFEM::addExternalForces() " << endl;
  // if(firstIter)     calcAndAssyLoadVector(1.0, 0.0);

  return;
}



int HBSplineFEM::factoriseSolveAndUpdate()
{
   //cout << "     HBSplineFEM: solving the matrix system ...  \n\n";
   char fct[] = "HBSplineFEM::factoriseSolveAndUpdate";
   //computerTime.go(fct);

   time_t tstart, tend;

   int ii, kk, ee;

  //cout << " rhsVec " << endl;        printVector(&(solverEigen->rhsVec(0)), totalDOF);
  //cout << " rhsVec " << endl;
  //for(int ii=fluidDOF-10;ii<totalDOF;ii++)
    //printf("%5d \t %12.8f \n", ii, solverEigen->rhsVec(ii));

  tstart = time(0);

  solverEigen->factoriseAndSolve();

   /*
   if(IterNum)
   {
     //cout << " factorising  " << endl;
     solver->factorise();
   }

   solver->solve();
   */
   
   tend = time(0);
   //printf("HBSplineFEM::factoriseSolveAndUpdate() took %8.4f second(s) \n ", difftime(tend, tstart) );

  //cout << " result " << endl;
  //for(int ii=fluidDOF;ii<totalDOF;ii++)
    //printf("%5d \t %12.8f \n", ii, solverEigen->soln(ii));

   
   double  *sln = &(solverEigen->soln[0]);

   for(ii=0;ii<velDOF;ii++)
     SolnData.var1[ii] += sln[ii];

   kk = velDOF;
   for(ii=0;ii<presDOF;ii++)
     SolnData.var2[ii] += sln[kk+ii];

   kk = velDOF + presDOF;
   for(ii=0;ii<IBDOF;ii++)
     SolnData.var3[ii] += sln[kk+ii];

   kk = velDOF + presDOF + IBDOF;
   for(ii=0;ii<solidDOF;ii++)
     SolnData.var4[ii] += sln[kk+ii];

   //computerTime.stopAndPrint(fct);

   //ctimFactSolvUpdt += computerTime.stop(fct);

   //cout << " Lagrange multipliers " << endl;        printVector(SolnData.var3);

    /*
    double  Fx=0.0, Fy=0.0;
    for(ii=0;ii<ImmersedBodyObjects[0]->nImmInt;ii++)
    {
      kk=ii*2;
      printf("%5d \t %20.16f \t %20.16f \n", ii, SolnData.var3(kk), SolnData.var3(kk+1));
      Fx += SolnData.var3(kk);
      Fy += SolnData.var3(kk+1);
    }
    
    printf("\n Forces = %20.16f \t %20.16f \n", Fx, Fy);
    */

   //cout << " result " << endl;        printVector(&(soln(0)), totalDOF);

   if(IBDOF > 111110)
   {
     cout << " result " << endl;//        printVector(&(soln(0)), totalDOF);
     for(int ii=fluidDOF-50;ii<totalDOF;ii++)
       printf("%5d \t %20.16f \n", ii, SolnData.var1(ii));
   }
   
   //solver->computeConditionNumber();

   //int ind=totalDOF/ndof, kk;
   //for(kk=0;kk<ind;kk++)	cout << kk << '\t' << soln(ndof*kk+2) << endl;

   //cout << "     HBSplineFEM: solving the matrix system ...DONE  \n\n";

   //computerTime.stopAndPrint(fct);

   return 0;
}



void HBSplineFEM::computeElementErrors(int index)
{
    int  ii, ee, count=0;
    node  *nd;

    totalError = 0.0;
    
    if(index < 4) // L2 or H1 norm based errors
    {
       for(ii=0;ii<activeElements.size();ii++)
       {
          ee = activeElements[ii];

          //cout << ee << '\t' << elems[ee]->getError() << endl;
          //totalError += ( elems[ee]->getError() * elems[ee]->getError() );
          totalError +=  elems[ee]->calcError(index);
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



void HBSplineFEM::computeConditionNumber()
{
  solverEigen->computeConditionNumber();

  return;
}





