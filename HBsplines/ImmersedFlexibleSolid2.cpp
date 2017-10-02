
#include "ImmersedFlexibleSolid.h"
#include "SolverMA41Eigen.h"
#include "MpapTime.h"
#include "TimeFunction.h"
#include "ImmersedIntegrationElement.h"
#include "../mySolvers/SolverPetsc.h"
#include "QuadratureUtil.h"


extern MpapTime mpapTime;
extern List<TimeFunction> timeFunction;



void ImmersedFlexibleSolid::initialise()
{
  SolnData.STAGGERED = STAGGERED;

  setSolver(1);
  setTimeParam();
  computeInitialAcceleration();

  return;
}




void ImmersedFlexibleSolid::setSolver(int slv, int *parm, bool cIO)
{
    //if(solver != NULL)
      //delete solver;
    //solver = NULL;

    int numProc=0;

    switch(slv)
    {
        case  1: // MA41 ..........................

            solver = (SolverEigen*) new SolverMA41Eigen;

            prepareMatrixPattern();

            //solver->printInfo();

            if(solver->initialise(0,0,totalDOF) != 0)
              return;

        break;

        case  4: // SolverEigen ..........................

            solver = (SolverEigen*) new SolverEigen;

            solver->setAlgorithmType(1);

            prepareMatrixPattern();

            if(solver->initialise(0,0,totalDOF) != 0)
              return;

            solver->printInfo();

        break;

        default: // invalid slv ...................

             cout << " this solver has not been implemented yet!\n\n";

        break;
    }

    solverOK = true;
    
    if(solver != NULL)
      solver->checkIO = cIO;

    //if(tis > 0)
      //setInitialConditions();

    for(int ii=0;ii<nElem;ii++)
    {
      elems[ii]->prepareElemData();
    }

    return;
}



void ImmersedFlexibleSolid::prepareMatrixPattern()
{
    //printf("\n     ImmersedFlexibleSolid::prepareMatrixPattern()  .... STARTED ...\n");

    int  r=0, c=0, r1=0, c1=0, count=0, count1=0, count2=0;
    int  ii=0, jj=0, iii=0, e=0, ee=0, ind=0, val1=0, val2=0;
    int  *tt, *tt1, *tt2, n1=0, n2=0, kk=0, e1=0, a=0, b=0, ll=0, pp=0, nnz=0;
    int  nRow=0, nCol=0, ind1=0, ind2=0, ndof_temp1=0;

    ndof = elems[0]->getNdofPerNode();
    ndof_temp1 = ndof;

    //totalDOF = nNode*ndof;     totalDOF += nElem_Constraint;

    node_map_new_to_old.resize(nNode, 0);
    node_map_old_to_new.resize(nNode, 0);

    dof_map_new_to_old.resize(totalDOF, 0);
    dof_map_old_to_new.resize(totalDOF, 0);

    kk=0;
    for(ii=0; ii<nNode; ii++)
    {
      node_map_new_to_old[ii] = ii;
      node_map_old_to_new[ii] = ii;

      for(jj=0; jj<ndof; jj++)
      {
        dof_map_new_to_old[kk] = kk;
        dof_map_old_to_new[kk] = kk;
        kk++;
      }
    }

    SolnData.node_map_new_to_old = node_map_new_to_old;
    SolnData.node_map_old_to_new = node_map_old_to_new;

    GeomData.node_map_new_to_old = node_map_new_to_old;
    GeomData.node_map_old_to_new = node_map_old_to_new;

    //cout << " nelem and  ndof " << nElem << '\t' << ndof << '\t' << npElem << endl;
    //cout << " ImmersedFlexibleSolid.... nNode and  ndof " << nNode << '\t' << ndof << endl;

    totalDOF = 0;
    for(ii=0;ii<nNode;ii++)
    {
      for(jj=0;jj<ndof;jj++)
      {
        //cout << ii << '\t' << jj << '\t' << NodeType[ii][jj] << endl;
        if(NodeType[ii][jj] == (int) -7777)
        {
          ID[ii][jj] = totalDOF++;

          assy4r.push_back(ii*ndof + jj);
        }
      }
    }

    // adjust assy4r array to account for the constraints

    ind = nNode*ndof_temp1;
    for(ii=0; ii<nElem_Constraint; ii++)
    {
      assy4r.push_back(ind+ii);
    }

    GeomData.assy4r = assy4r;

    // adjust ID array to account for the constraints

    for(ee=0; ee<nElem; ee++)
    {
      ii = elems[ee]->nodeNums[0];

      pp = elems[ee]->getElmTypeNameNum();

      if( pp == 33 )
      {
        jj = ndof_temp1 + 1;
        ID[ii][jj] = totalDOF++;
      }

      if( pp == 34 )
      {
        jj = ndof_temp1 + 0;
        ID[ii][jj] = totalDOF++;
      }
    }

    //cout << " totalDOF " << totalDOF << endl;
    //assy4r.resize(totalDOF);


    for(ee=0; ee<nElem; ee++)
    {
      //cout << " ee = " << ee << endl;
      npElem     = elems[ee]->getNodesPerElement();
      ndof_temp1 = elems[ee]->getNdofPerNode();
      ind        = ndof_temp1*npElem;

      //printVector(IEN[ee]);

      pp = elems[ee]->getElmTypeNameNum();

      if( pp == 33 ) // constaint element with rigid X-axis
      {
        LM[ee].resize(2);  // to account for the Lagrange multiplier

        kk = IEN[ee][0];

        LM[ee][0] = ID[kk][1];
        LM[ee][1] = ID[kk][3];
      }
      else if( pp == 34 ) // constaint element with rigid Y-axis
      {
        LM[ee].resize(2);  // to account for the Lagrange multiplier

        kk = IEN[ee][0];

        LM[ee][0] = ID[kk][0];
        LM[ee][1] = ID[kk][2];
      }
      else // standard elements
      {
        LM[ee].resize(ind);
        for(ii=0;ii<npElem;ii++)
        {
          ind = ndof_temp1*ii;

          kk = IEN[ee][ii];

          for(jj=0;jj<ndof_temp1;jj++)
          {
            //cout << ee << '\t' << ii << '\t' << jj << '\t' << ind << '\t' << ID[kk][jj] << endl;
            //cout << ee << '\t' << ii << '\t' << jj << '\t' << ind << '\t' << LM[ee][ind+jj] << '\t' << ID[IEN[ee][ii]][jj] << endl;
            LM[ee][ind+jj] = ID[kk][jj];
            //cout << " IIIIIIIII " << endl;
          }
        }
      }
      //cout << " ee = " << ee << endl;    
    }

    //cout << " totalDOF " << totalDOF << endl;
    //cout << " nElem " << nElem << endl;
    for(ee=0; ee<nElem; ee++)
    {
      elems[ee]->forAssyVec = LM[ee];
    }

    pp=false;
    //pp=true;
    if(pp)
    {
       printf("   IEN array \n\n");
       for(ee=0; ee<nElem; ee++)
       {
          npElem = elems[ee]->getNodesPerElement();
          for(ii=0; ii<npElem; ii++)
            cout << '\t' << IEN[ee][ii];
          cout << endl;
       }
       printf("\n\n\n");

       printf("   ID array \n\n");
       for(ii=0; ii<nNode; ii++)
       {
          for(jj=0; jj<ID[ii].size(); jj++)
            cout << '\t' << ID[ii][jj];
          cout << endl;
       }
       printf("\n\n\n");

       printf("   LM array \n\n");
       for(ee=0; ee<nElem; ee++)
       {
          for(jj=0; jj<LM[ee].size(); jj++)
            cout << '\t' << LM[ee][jj];
          cout << endl;
       }
       printf("\n\n\n");

       printf("  assy4r array \n\n");
       for(ii=0;ii<totalDOF;ii++)
       {
          cout << ii << '\t' << assy4r[ii] << endl;
       }
       printf("\n\n\n");
    }

    //printf("\n element DOF values initialised \n\n");
    //printf("\n Preparing matrix pattern \n\n");

    vector<int>::const_iterator location;
    set<int>::iterator it;

    forAssyMat.resize(totalDOF);
    
    for(ee=0;ee<nElem;ee++)
    {
       tt = &(LM[ee][0]);
       nsize = LM[ee].size();

       for(ii=0;ii<nsize;ii++)
       {
          count1 = tt[ii];

          if(tt[ii] != -1)
          {
            for(jj=0;jj<nsize;jj++)
            {
              if(tt[jj] != -1)
                forAssyMat[count1].push_back(tt[jj]);
            }
          }
       }
    }

    bool pp1=false;
    //pp1=true;
    if(pp1)
    {
       printf("   Number of non-zeros = %5d \n\n", nnz);
       printf("   dof to dof connectivity ...:  \n\n");
       for(ii=0;ii<totalDOF;ii++)
       {
          cout << " dof # " << ii << " : ";
          for(jj=0;jj<forAssyMat[ii].size();jj++)
            cout << '\t' << forAssyMat[ii][jj];
          cout << endl;
       }
       printf("\n\n\n");
    }

    //printf("\n Preparing matrix pattern DONE \n\n");

    VectorXd  nnzVec(totalDOF);

    nnz = 0;
    for(ii=0;ii<totalDOF;ii++)
    {
      forAssyMat[ii].push_back(ii);
      findUnique(forAssyMat[ii]);
      nnzVec[ii] = forAssyMat[ii].size();
      nnz += nnzVec[ii];
    }
    //cout << " nnz " << nnz << endl;

    nRow = nCol = totalDOF;

    //cout << " AAAAAAAAAA " << endl;
    solver->mtx.setZero();

    solver->mtx.resize(nRow, nCol);
    solver->mtx.reserve(nnz);
    solver->mtx.reserve(nnzVec);

    for(ii=0;ii<totalDOF;ii++)
    {
      for(jj=0;jj<forAssyMat[ii].size();jj++)
      {
        //cout << ii << '\t' << forAssyMat[ii][jj] << endl;
        solver->mtx.coeffRef(ii, forAssyMat[ii][jj]) = 0.0;
      }
    }

    solver->mtx.makeCompressed();

    solver->currentStatus = PATTERN_OK;

    if(!STAGGERED)
    {
      forAssyCoupledHorz.resize(nNode*ndof_temp1);
      forAssyCoupledVert.resize(nNode*DIM);

      for(ii=0;ii<nNode;ii++)
      {
        ind1 = ii*ndof_temp1;
        ind2 = ii*DIM;

        for(jj=0;jj<DIM;jj++)
        {
          forAssyCoupledHorz[ind1+jj].push_back(ind2+jj);
          //forAssyCoupledVert[ind2+jj].push_back(ind1+jj);
        }
      }
      //for(ii=0;ii<totalDOF;ii++)
        //printVector(forAssyCoupledHorz[ii]);
    }

    //printf("\n     ImmersedFlexibleSolid::prepareMatrixPattern()  .... FINISHED ...\n\n");

    return;
}




void ImmersedFlexibleSolid::solveTimeStep()
{
    printf("\n Solving Immersed Flexible Solid \n");
    //printf("\n External force norm = %12.6E \n", forceCur.norm());

    int  ii=0, ee=0;

    //printVector(SolnData.forceCur);

    for(ii=0;ii<10;ii++)
    {
      calcStiffnessAndResidual(1, 1, 1);

      printf("\t %5d \t %12.6E\n", ii, rNorm);
      //printf("\t %5d \t %5d \t %12.6E\n", ii, firstIter, rNorm);

      if(converged())
        break;

      factoriseSolveAndUpdate();

      updateIterStep();
    }

    printf("\n Solving Immersed Flexible Solid ..... DONE  \n\n");

    return;
}






int ImmersedFlexibleSolid::calcStiffnessAndResidual(int printRes, bool zeroMtx, bool zeroRes)
{
    //cout << "     ImmersedFlexibleSolid: generating coefficient Matrices ...\n\n";

    if(solver == NULL)
    {
      COUT << "You need to select a solver first!\n\n";
      return 1;
    }

    solver->zeroMtx();

    if(firstIter)
      rNorm = -1.0;

    MatrixXd  Klocal;
    VectorXd  Flocal;

    for(int ee=0;ee<nElem;ee++)  // loop over all the elements
    {
      //cout << "       elem... : " << (ee+1) << endl;

      elems[ee]->calcStiffnessAndResidual(Klocal, Flocal);

      //cout << " MMMMMMMMMMM " << endl;
      //elems[ee]->assembleElementMatrixAndVector(0, solver->mtx, &(solver->rhsVec(0)));

      //elems[ee]->assembleElementMatrix(0, solver->mtx);
      //cout << " MMMMMMMMMMM " << endl;
      //elems[ee]->assembleElementVector(false, false, &(solver->rhsVec(0)), &(SolnData.reac(0)), 0, 0);

      solver->assembleMatrixAndVector(0, 0, elems[ee]->forAssyVec, Klocal, Flocal);
    }

    //cout << " solver->rhsVec " << endl;        printVector(solver->rhsVec);

    //printf("\n rhsVec norm = %12.6E \n", solver->rhsVec.norm());

    //applyBoundaryConditions(0, 0, solver->mtx, &(solver->rhsVec(0)));

    applyExternalForces();

    // rhs due to external forces
    //rhsVec += rhsVec2;

    //solver->rhsVec(totalDOF-2) += 0.5;

    //cout << " rhsVec " << endl;        printVector(&(rhsVec[0]), totalDOF);

    //printf("\n rhsVec norm = %12.6E \n", solver->rhsVec.norm());

    firstIter = false;
    SolnData.firstIter = firstIter;
    rNormPrev = rNorm;
    rNorm     = solver->rhsVec.norm();

    if (printRes > 1) { COUT << "ImmersedFlexibleSolid"; printf("  %11.4e\n",rNorm);}

    solver->currentStatus = ASSEMBLY_OK;

    return 0;
}




int ImmersedFlexibleSolid::applyBoundaryConditions(int start1, int start2, SparseMatrixXd& globalK, double* rhs)
{   
    int  ii=0, jj=0, nn=0, dof=0, aa=0, ind=0;
    double  specVal=0.0, PENALTY=1.0e8;
    double  af = SolnData.td(2);

    for(aa=0;aa<DirichletBCs.size();aa++)
    {
      nn  = (int) (DirichletBCs[aa][0] - 1);
      dof = (int) (DirichletBCs[aa][1] - 1);
      specVal = DirichletBCs[aa][2];

      ind = nn*ndof+dof;
      //specVal  -=  SolnData.var1DotCur[ind];
      specVal  -=  SolnData.var1Cur[ind];

      //cout << start << '\t' << nn << '\t' << ind << '\t' << specVal << endl;

      ind += start1;
      globalK.coeffRef(ind, ind) += af*PENALTY;
      rhs[ind]   += (PENALTY*specVal);
    }

    return 0;
}



int  ImmersedFlexibleSolid::applyExternalForces()
{
    //printVector(SolnData.forceCur);

    for(int ii=0;ii<totalDOF;ii++)
    {
      solver->rhsVec[ii] += SolnData.forceCur[assy4r[ii]];
      //solver->rhsVec[ii] += fluidAcceCur[assy4r[ii]];
    }

  /*
  ///////////////////
  // stiffness contribution from PENALTY

  double  dvol, detJ, fact, dt=mpapTime.dt;
  double  velFact = SolnData.td[6];

  int  nlb = ImmIntgElems[0]->pointNums.size();

  int  nGauss = 5;

  VectorXd  Nb;    Nb.resize(nlb);
  myPoint  geom, param;

  vector<double>  gausspoints, gaussweights;

  getGaussPoints1D(nGauss, gausspoints, gaussweights);

  int  aa, ii, jj, ind, dd, gp;

    ImmersedIntegrationElement  *lme;
    myPoly*  poly;


    for(aa=0; aa<nImmInt; aa++)
    {
      lme  = ImmIntgElems[aa];
      poly = ImmersedFaces[aa];

      //cout << " aa = " << aa << '\t' << lme->isActive() << endl;
      //if( lme->isActive() )
      //{
          for(gp=0; gp<nGauss; gp++)
          {
              param[0] = gausspoints[gp];

              poly->computeBasisFunctions(param, geom, Nb, detJ);

              dvol  = gaussweights[gp] * detJ;

              //printf(" %12.6f,  %12.6f, %12.6f \n", Nb[0], Nb[1], dvol);

              dvol *= PENALTY;
              for(ii=0; ii<nlb; ii++)
              {
                ind = ImmIntgElems[aa]->pointNums[ii]*2;

                for(dd=0; dd<2; dd++)
                {
                  solver->rhsVec[ind+dd]  +=  (SolnData.var1DotPrev[ind+dd]-SolnData.var1DotCur[ind+dd])*dvol;
                  solver->mtx.coeffRef(ind+dd, ind+dd) += velFact*dvol;
                }
              }
          }
      //}
    }
    */

  return 0;
}



void ImmersedFlexibleSolid::calcForceVector()
{
    //for(int ee=0;ee<nElem;ee++)
      //elems[ee]->calcExternalForces();

    return;
}




int ImmersedFlexibleSolid::factoriseSolveAndUpdate()
{
    //cout << " residue_new " << endl;        printVector(&(solver->rhsVec[0]), totalDOF);

    solver->factoriseAndSolve();

    soln.setZero();
    // update solution vector
    for(int kk=0;kk<totalDOF;kk++)
      soln[assy4r[kk]] = solver->soln[kk];

    SolnData.var1 += soln;
    //SolnData.var1Dot += soln;

    //printVector(SolnData.var1);

    //printf("\n\n\n");

    //cout << " result " << endl;        printVector(&(SolnData.var1[0]), totalDOF);

    //ctimFactSolvUpdt += computerTime.stop(fct);

    return 0;
}



void ImmersedFlexibleSolid::calcCouplingMatrices()
{
    //////////////////////////////////////////
    // off-diagonal matrices
    //////////////////////////////////////////

    int  aa=0, nlb=0, ii=0, jj=0, kk=0, nlbS=0, nlbL=0;
    int  TI=0, TIp1=0, TIp2=0, TJ=0, TJp1=0;
    MatrixXd  Ktemp, Ktemp2;

    int  ind1 = nNode*ndof;
    int  ind2 = nNode*DIM;

    Khorz.resize(ind1, ind2);
    Khorz.setZero();

    Kvert.resize(ind2, ind1);
    Kvert.setZero();

    for(aa=0;aa<ImmIntgElems.size();aa++)
    {
      nlbL = ImmIntgElems[aa]->pointNums.size();
      nlbS = elems[aa]->nodeNums.size();

      //cout << " aa = " << aa << endl;
      ImmIntgElems[aa]->computeKhorzKvertFlexible(0, 0, Ktemp, Ktemp2);
      //printMatrix(Ktemp);
      //printf("\n\n");
      //printMatrix(Ktemp2);
      //printf("\n\n");

      for(ii=0; ii<nlbS; ii++)
      {
        TI   = 3*elems[aa]->nodeNums[ii];
        TIp1 = TI+1;
        TIp2 = TI+2;
        
        ind1 = 3*ii;

        for(jj=0; jj<nlbL; jj++)
        {
          TJ   = ImmIntgElems[aa]->pointNums[jj] * DIM;
          TJp1 = TJ+1;
          
          ind2 = jj*DIM;

          Khorz(TI, TJ)      +=  Ktemp(ind1, ind2);
          Khorz(TI, TJ+1)    +=  Ktemp(ind1, ind2+1);

          Khorz(TI+1, TJ)    +=  Ktemp(ind1+1, ind2);
          Khorz(TI+1, TJ+1)  +=  Ktemp(ind1+1, ind2+1);

          Kvert(TJ,   TI)    +=  Ktemp2(ind2, ind1);
          Kvert(TJ+1, TI)    +=  Ktemp2(ind2+1, ind1);

          Kvert(TJ,   TI+1)  +=  Ktemp2(ind2, ind1+1);
          Kvert(TJ+1, TI+1)  +=  Ktemp2(ind2+1, ind1+1);
        }
      }
    }

    Khorz *= -SolnData.td[2];
    Kvert *= -SolnData.td[2];

  //printMatrix(Khorz);
}





int ImmersedFlexibleSolid::assembleGlobalMatrixAndVector(int start1, int start2, SparseMatrixXd& mtx, double* rhs)
{
    if(totalDOF <= 0)
      return 1;

    int ee=0, ii=0, jj=0, k1=0, k2=0, r=0, c=0;

    for(ee=0;ee<nElem;ee++)  // loop over all the elements
    {
      //cout << "       elem... : " << (ee+1) << endl;

      //elems[ee]->resetMatrixAndVector();

      elems[ee]->calcStiffnessAndResidual();

      //cout << " MMMMMMMMMMM " << endl;
      elems[ee]->assembleElementMatrixAndVector(start2, mtx, rhs);
    }

    //cout << " aaaaaaaaaaaaaaa " << endl;
    for(ii=0;ii<totalDOF;ii++)
      rhs[start2+ii] += SolnData.forceCur[assy4r[ii]];

    applyBoundaryConditions(start1, start2, mtx, rhs);
  
    //calcCouplingMatrices();

//
    double af = SolnData.td[2];

    for(ii=0;ii<totalDOF;ii++)
    {
      r  = start2 + ii;
      k1 = assy4r[ii];

      for(jj=0;jj<forAssyCoupledHorz[ii].size();jj++)
      {
        c = start1 + forAssyCoupledHorz[ii][jj];

        mtx.coeffRef(r, c) += -af;
        mtx.coeffRef(c, r) += -af;
      }
    }
//

/*
    for(ii=0;ii<totalDOF;ii++)
    {
      r  = start2 + ii;
      k1 = assy4r[ii];

      for(jj=0;jj<forAssyCoupledHorz[ii].size();jj++)
      {
        k2 = forAssyCoupledHorz[ii][jj];
        c = start1 + k2;

        mtx.coeffRef(r, c) += Khorz(k1, k2);
        mtx.coeffRef(c, r) += Kvert(k2, k1);
      }
    }
*/

  return 0;
}


/*
int ImmersedFlexibleSolid::AppendSolidMatrixPatternCutFEM(int start1, int start2, vector<vector<int> >& forAssyGlobal)
{
  for(ii=0; ii<forAssyMat.size(); ii++)
  {
    for(jj=0; jj<forAssyMat[ii].size(); jj++)
    {
      forAssyGlobal[ii].push_back(forAssyMat[ii][jj]);
    }
  }

  return 0;
}
*/


int ImmersedFlexibleSolid::assembleGlobalMatrixAndVectorCutFEM(int start1, int start2, SolverPetsc* solverTemp)
{
    int ee=0, ii=0, jj=0, ind=0;

    MatrixXd  Klocal;
    VectorXd  Flocal;

    for(ee=0; ee<nElem; ee++)  // loop over all the elements
    {
      //cout << "       elem... : " << (ee+1) << endl;

      elems[ee]->calcStiffnessAndResidual(Klocal, Flocal);

      solverTemp->assembleMatrixAndVector(start1, start2, elems[ee]->forAssyVec, elems[ee]->forAssyVec, Klocal, Flocal);
    }

    //cout << " solver->rhsVec " << endl;        printVector(solver->rhsVec);

    //printf("\n rhsVec norm = %12.6E \n", solver->rhsVec.norm());

    applyBoundaryConditions(start1, start2, solverTemp->mtx, solverTemp->rhsVec);

    //applyExternalForces();

    //cout << " aaaaaaaaaaaaaaa " << endl;
    for(ii=0;ii<totalDOF;ii++)
    {
      ind = start2+ii;
      VecSetValue(solverTemp->rhsVec, ind, SolnData.forceCur[assy4r[ii]], ADD_VALUES);
    }

    return 0;
}




int ImmersedFlexibleSolid::applyBoundaryConditions(int start1, int start2, Mat mtxTemp, Vec rhsTemp)
{   
    int  ii=0, jj=0, nn=0, dof=0, aa=0, ind=0;
    double  specVal=0.0, PENALTY=1.0e8, stiff=0.0, res=0.0;

    double  af = SolnData.td(2);
    double  velFact = SolnData.td(10);

    for(aa=0;aa<DirichletBCs.size();aa++)
    {
      nn  = (int) (DirichletBCs[aa][0] - 1);
      dof = (int) (DirichletBCs[aa][1] - 1);
      specVal = DirichletBCs[aa][2];

      ind = nn*ndof+dof;
      //specVal  -=  SolnData.var1DotCur[ind];
      specVal  -=  SolnData.var1Cur[ind];

      //cout << start << '\t' << nn << '\t' << ind << '\t' << specVal << endl;

      ind += start1;

      stiff = af*PENALTY/velFact;
      res   = specVal*PENALTY;

      MatSetValue(mtxTemp, ind, ind, stiff, ADD_VALUES);
      VecSetValue(rhsTemp, ind, res, ADD_VALUES);
    }

    return 0;
}



