
#include "ImmersedFlexibleSolid.h"

#include "SolverMA41Eigen.h"
#include "MpapTime.h"
#include "TimeFunction.h"
#include "ImmersedIntegrationElement.h"


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

    char fct[] = "ImmersedFlexibleSolid::setSolver";
    
    int numProc;
    
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
  printf("\n     ImmersedFlexibleSolid::prepareMatrixPattern()  .... STARTED ...\n");

  char fct[] = "ImmersedFlexibleSolid::prepareMatrixPattern";

    int  r, c, r1, c1, count=0, count1=0, count2=0, ii, jj, iii, e, ee, ind;
    int *tt, *tt1, *tt2,  val1, val2, n1, n2, kk, e1, a, b, ll, pp, nnz;
    int nRow, nCol, ind1, ind2;


    totalDOF = nNode*ndof;

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

    cout << " nelem and  ndof " << nElem << '\t' << ndof << '\t' << npElem << endl;

    cout << " ImmersedFlexibleSolid.... nNode and  ndof " << nNode << '\t' << ndof << endl;


    totalDOF = 0;
    for(ii=0;ii<nNode;ii++)
    {
      for(jj=0;jj<ndof;jj++)
      {
        //cout << ii << '\t' << jj << '\t' << NodeType[ii][jj] << endl;
        if(NodeType[ii][jj] == (int) -7777)
        {
          ID[ii][jj] = totalDOF++;
        }
      }
    }

    cout << " totalDOF " << totalDOF << endl;
    assy4r.resize(totalDOF);

    ind = ndof*npElem;

    IEN.resize(nElem);
    LM.resize(nElem);
    for(ee=0; ee<nElem; ee++)
    {
      IEN[ee] = elems[ee]->nodeNums;
      LM[ee].resize(ind);
    }

    cout << ndof << '\t' << npElem << '\t' << ind << '\t' << nElem << endl;

    for(ee=0;ee<nElem;ee++)
    {
      for(ii=0;ii<npElem;ii++)
      {
        ind = ndof*ii;

        kk = IEN[ee][ii];

        for(jj=0;jj<ndof;jj++)
        {
          //cout << ee << '\t' << ii << '\t' << jj << '\t' << ind << '\t' << ID[kk][jj] << endl;
          //cout << ee << '\t' << ii << '\t' << jj << '\t' << ind << '\t' << LM[ee][ind+jj] << '\t' << ID[IEN[ee][ii]][jj] << endl;
          LM[ee][ind+jj] = ID[kk][jj];
          //cout << " IIIIIIIII " << endl;
        }
      }
    }

    cout << " totalDOF " << totalDOF << endl;
    count = 0;
    for(ii=0;ii<nNode;ii++)
    {
      for(jj=0;jj<ndof;jj++)
      {
        //cout << ii << '\t' << jj << '\t' << ID[ii][jj] << endl;
        if( ID[ii][jj] != -1)
          assy4r[count++] = ii*ndof + jj;
      }
    }
    cout << " totalDOF " << totalDOF << endl;
    cout << " nElem " << nElem << endl;
    for(ii=0;ii<nElem;ii++)
    {
      elems[ii]->forAssyVec = LM[ii];
    }

    cout << " totalDOF " << totalDOF << endl;
    pp=false;
    //pp=true;
    if(pp)
    {
       printf("   IEN array \n\n");
       for(ii=0;ii<nElem;ii++)
       {
          for(jj=0;jj<npElem;jj++)
            cout << '\t' << IEN[ii][jj];
          cout << endl;
       }
       printf("\n\n\n");

       printf("   ID array \n\n");
       for(ii=0;ii<nNode;ii++)
       {
          for(jj=0;jj<ndof;jj++)
            cout << '\t' << ID[ii][jj];
          cout << endl;
       }
       printf("\n\n\n");

       printf("   LM array \n\n");
       for(ii=0;ii<nElem;ii++)
       {
          for(jj=0;jj<nsize;jj++)
            cout << '\t' << LM[ii][jj];
          cout << endl;
       }
       printf("\n\n\n");

       printf("  assy4r array \n\n");
       for(ii=0;ii<totalDOF;ii++)
       {
          cout << assy4r[ii] << endl;
       }
       printf("\n\n\n");
    }

    printf("\n element DOF values initialised \n\n");
    printf("\n Preparing matrix pattern \n\n");

    vector<int>::const_iterator location;
    set<int>::iterator it;

    forAssyMat.resize(totalDOF);
    
    for(ee=0;ee<nElem;ee++)
    {
       tt = &(LM[ee][0]);

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

    printf("\n Preparing matrix pattern DONE \n\n");

    VectorXd  nnzVec(totalDOF);

    nnz = 0;
    for(ii=0;ii<totalDOF;ii++)
    {
      forAssyMat[ii].push_back(ii);
      findUnique(forAssyMat[ii]);
      nnzVec[ii] = forAssyMat[ii].size();
      nnz += nnzVec[ii];
    }
    cout << " nnz " << nnz << endl;

    nRow = nCol = totalDOF;

    cout << " AAAAAAAAAA " << endl;
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
      forAssyCoupledHorz.resize(nNode*ndof);
      forAssyCoupledVert.resize(nNode*DIM);

      for(ii=0;ii<nNode;ii++)
      {
        ind1 = ii*ndof;
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

    printf("\n     ImmersedFlexibleSolid::prepareMatrixPattern()  .... FINISHED ...\n\n");

    return;
}




void ImmersedFlexibleSolid::SolveTimeStep()
{
  printf("\n Solving Immersed Flexible Solid \n");
  //printf("\n External force norm = %12.6E \n", forceCur.norm());

  int  ii, ee;

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
      //elems[ee]->AssembleElementMatrixAndVector(0, solver->mtx, &(solver->rhsVec(0)));

      //elems[ee]->AssembleElementMatrix(0, solver->mtx);
      //cout << " MMMMMMMMMMM " << endl;
      //elems[ee]->AssembleElementVector(false, false, &(solver->rhsVec(0)), &(SolnData.reac(0)), 0, 0);

      solver->AssembleMatrixAndVector(0, 0, elems[ee]->forAssyVec, Klocal, Flocal);
  }

  //cout << " solver->rhsVec " << endl;        printVector(solver->rhsVec);

  //printf("\n rhsVec norm = %12.6E \n", solver->rhsVec.norm());

  //applyBoundaryConditions(0, solver->mtx, &(solver->rhsVec(0)));

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




void ImmersedFlexibleSolid::applyBoundaryConditions(int start, SparseMatrixXd& globalK, double* rhs)
{   
  int ii, jj, nn, dof, aa, ind;
  double specVal, PENALTY=1.0e8, af;

  af = SolnData.td(2);

  for(aa=0;aa<DirichletBCs.size();aa++)
  {
      nn  = (int) (DirichletBCs[aa][0] - 1);
      dof = (int) (DirichletBCs[aa][1] - 1);
      specVal = DirichletBCs[aa][2];

      ind = nn*ndof+dof;
      //specVal  -=  SolnData.var1DotCur[ind];
      specVal  -=  SolnData.var1Cur[ind];

      //cout << start << '\t' << nn << '\t' << ind << '\t' << specVal << endl;

      ind += start;
      globalK.coeffRef(ind, ind) += af*PENALTY;
      rhs[ind]   += (PENALTY*specVal);
  }

  return;
}



void ImmersedFlexibleSolid::applyExternalForces()
{
  //printVector(SolnData.forceCur);

  for(int ii=0;ii<totalDOF;ii++)
    solver->rhsVec[ii] += SolnData.forceCur[assy4r[ii]];

  return;
}



void ImmersedFlexibleSolid::calcForceVector()
{
  //for(int ee=0;ee<nElem;ee++)
    //elems[ee]->calcExternalForces();

  return;
}


int ImmersedFlexibleSolid::factoriseSolveAndUpdate()
{
  char fct[] = "ImmersedFlexibleSolid::factoriseSolveAndUpdate";
  //computerTime.go(fct);

  //cout << " residue_new " << endl;        printVector(&(solver->rhsVec[0]), totalDOF);

  solver->factoriseAndSolve();

  soln.setZero();
  // update solution vector
  for(int kk=0;kk<totalDOF;kk++)
    soln[assy4r[kk]] = solver->soln[kk];

  SolnData.var1 += soln;
  //SolnData.var1Dot += soln;

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

    int  aa, nlb, ii, jj, kk, ind1, ind2, nlbS, nlbL;
    int  TI, TIp1, TIp2, TJ, TJp1;
    MatrixXd  Ktemp, Ktemp2;

    ind1 = nNode*ndof;
    ind2 = nNode*DIM;

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





int ImmersedFlexibleSolid::AssembleGlobalMatrixAndVector(int start1, int start2, SparseMatrixXd& mtx, double* rhs)
{
  if(totalDOF <= 0)
    return 1;

  int ee, ii, jj, k1, k2, r, c;

  for(ee=0;ee<nElem;ee++)  // loop over all the elements
  {
    //cout << "       elem... : " << (ee+1) << endl;

    //elems[ee]->resetMatrixAndVector();

    elems[ee]->calcStiffnessAndResidual();

    //cout << " MMMMMMMMMMM " << endl;
    elems[ee]->AssembleElementMatrixAndVector(start2, mtx, rhs);
  }

  //cout << " aaaaaaaaaaaaaaa " << endl;
  for(ii=0;ii<totalDOF;ii++)
    rhs[start2+ii] += SolnData.forceCur[assy4r[ii]];

  applyBoundaryConditions(start2, mtx, rhs);
  
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



int ImmersedFlexibleSolid::AssembleGlobalMatrixAndVectorCutFEM(int start1, int start2, SparseMatrixXd& mtx, double* rhs)
{
  int ee, ii, jj, kk, r, c;
  for(ee=0;ee<nElem;ee++)  // loop over all the elements
  {
    //cout << "       elem... : " << (ee+1) << endl;

    //elems[ee]->resetMatrixAndVector();
    //elems[ee]->calcStiffnessAndResidual();

    //cout << " MMMMMMMMMMM " << endl;
    elems[ee]->AssembleElementMatrixAndVector(start2, mtx, rhs);
  }

  //cout << " aaaaaaaaaaaaaaa " << endl;
  for(ii=0;ii<totalDOF;ii++)
    rhs[start2+ii] += SolnData.forceCur[assy4r[ii]];

  applyBoundaryConditions(start2, mtx, rhs);
  
  double af = SolnData.td[2];

  for(ii=0;ii<totalDOF;ii++)
  {
    r  = start2 + ii;
    kk = assy4r[ii];

    for(jj=0;jj<forAssyCoupledHorz[ii].size();jj++)
    {
      c = start1 + forAssyCoupledHorz[ii][jj];

      mtx.coeffRef(r, c) += -af;
      mtx.coeffRef(c, r) += -af;
    }
  }

  return 0;
}



