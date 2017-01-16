
#include <iostream>
#include <bits/algorithmfwd.h>

#include "SolverMA41Eigen.h"
#include "FunctionsSolver.h"
#include "FunctionsProgram.h"
#include "SolverWorkSpace.h"
#include "SolverTime.h"
#include "ComputerTime.h"


extern SolverWorkSpace solverWorkSpace;
extern SolverTime      solverTime;
extern ComputerTime    computerTime;


using namespace std;


SolverMA41Eigen::SolverMA41Eigen()
{
  //comprMtxFlg = true;

  ROWSCA = new double [10];
  COLSCA = new double [10];

  IS     = NULL;
  
  numSCA = 10;

  return;
}




SolverMA41Eigen::~SolverMA41Eigen()
{
  free();
	
  if (ROWSCA != NULL) delete [] ROWSCA;
  if (COLSCA != NULL) delete [] COLSCA;
  if (IS     != NULL) delete [] IS;

  return;
}



int SolverMA41Eigen::initialise(int numProc, int matrixType, int nr)
{
  nRow = nCol = nr;

  soln.resize(nRow);
  soln.setZero();

  rhsVec   = soln;
  solnPrev = soln;

  //if (currentStatus != PATTERN_OK)
    //{ prgWarning(1,"SolverMA41Eigen::initialise","prepare matrix pattern first!"); return 1; }

  int  *tmp, ii, k, MAXS, JOB, count, i, m;

  double *S;

  N  = nRow;
  NE = mtx.nonZeros();

  // initialise controls

  ma41id_(CNTL,ICNTL,KEEP);

  CNTL [0] = 0.01;//.01;  // pivoting
//  ICNTL[0] = 1;
  //ICNTL[7] = 6;    // scaling
  //ICNTL[9] = 10;   // max iterations

//    for(int ii=0;ii<10;ii++)       cout << " ICNTL[" << ii << "] .... : " << ICNTL[ii] << endl;  cout << endl;

//  ICNTL[4] = 3;
  
  // analyse matrix pattern

  MAXIS = 2*NE + 11*N + 1;
  
  IS    = new int [MAXIS];
  
  JOB   = 1;

  MAXS  = solverWorkSpace.dimDbl;
  S     = solverWorkSpace.dbl;
  
  row.resize(NE);
  col.resize(NE);
  array.resize(NE);
  
  cout << nRow << '\t' << NE << endl;
  
  count = 0;
  for(k=0; k<mtx.outerSize(); ++k)
  {
    for(SparseMatrixXd::InnerIterator it(mtx,k); it; ++it)
    {
      row[count]   = it.row() + 1;
      col[count]   = it.col() + 1;
      array[count] = it.value();
      count++;
    }
  }

  ma41ad_(&JOB, &N, &NE, &row[0], &col[0], &array[0], &rhsVec[0], COLSCA,ROWSCA,KEEP,
	  IS,&MAXIS,S,&MAXS, CNTL,ICNTL,INFO,RINFO);

  // allocate memory for IS if required

//  m = INFO[6] + 10 + roundToInt(0.2 * (double)NE);

  m = roundToInt(5.0 * (double)INFO[6]);
  
//    cout << " MAXIS  && MAXS && m .... : " << MAXIS << '\t' << MAXS << '\t' << m << endl;
    
//    for(int ii=0;ii<10;ii++)       cout << " INFO[" << ii << "] .... : " << INFO[ii] << endl;
 
  if(MAXIS < m) 
  {
    tmp = new int [m];
    for (i=0; i<MAXIS; i++)
      tmp[i] = IS[i];

    delete [] IS; 
    IS = tmp;
    MAXIS = m;
  }
  
  // allocate memory for S if required
  
  m = roundToInt(10.0 * (double)INFO[7]);
  
//  cout << " m .... : " << m << endl;
//  cout << " KEEP[12] .... : " << KEEP[12] << endl;
  
  solverWorkSpace.expand(m);

  currentStatus = INIT_OK;

  return 0;
}








int SolverMA41Eigen::factorise()
{
  char fct[] = "SolverMA41Eigen::factorise";

  if(currentStatus != ASSEMBLY_OK) { prgWarning(1,fct,"assemble matrix first!"); return 1; }

  int  JOB  = 2, MAXS = solverWorkSpace.dimDbl;

  double *S = solverWorkSpace.dbl;

  if (checkIO)
  {
    // search for "nan" entries in matrix coefficients
    //if (prgNAN(mtx.x.x,NE)) prgError(1,fct,"nan matrix coefficient!");
  }

  computerTime.go(fct);

  ma41ad_(&JOB,&N,&NE,&row[0], &col[0], &array[0], &rhsVec[0],COLSCA,ROWSCA,KEEP,IS,&MAXIS,S,&MAXS,
               CNTL,ICNTL,INFO,RINFO);
 
  currentStatus = FACTORISE_OK;
  
  solverTime.total     -= solverTime.factorise;
  solverTime.factorise += computerTime.stop(fct);
  solverTime.total     += solverTime.factorise;

  return INFO[0];
}






int SolverMA41Eigen::solve()
{
  char fct[] = "SolverMA41Eigen::solve";

  if (currentStatus != FACTORISE_OK) { prgWarning(1,fct,"factorise matrix first!"); return -1; }

  int    JOB  = 3, MAXS = solverWorkSpace.dimDbl;

  double *S = solverWorkSpace.dbl;

  if(checkIO)
  {
    // search for "nan" entries in RHS
    //if (prgNAN(RHS,N)) prgError(1,fct,"nan rhs entry!");
  }

  computerTime.go(fct);

  ma41ad_(&JOB,&N,&NE,&row[0], &col[0], &array[0], &rhsVec[0],COLSCA,ROWSCA,KEEP,IS,&MAXIS,S,&MAXS,
             CNTL,ICNTL,INFO,RINFO);

  soln = rhsVec;

  solverTime.total -= solverTime.solve;
  solverTime.solve += computerTime.stop(fct);
  solverTime.total += solverTime.solve;

  if(INFO[0] != 0)
    return 2;
  
  if(checkIO)
  {
    // search for "nan" entries in solution vector
    //if (prgNAN(RHS,N)) prgError(1,fct,"nan entry in solution vector!");
  }

  return 0;
}






int SolverMA41Eigen::factoriseAndSolve()
{
  char fct[] = "SolverMA41Eigen::factoriseAndSolve";

  if (currentStatus != ASSEMBLY_OK) { prgWarning(1,fct,"assemble matrix first!"); return -1; }

  int    JOB  = 5, MAXS = solverWorkSpace.dimDbl, k, count;

  double *S = solverWorkSpace.dbl;

  if(checkIO)
  {
    // search for "nan" entries in matrix coefficients
    //if (prgNAN(mtx.x.x,NE)) prgError(1,fct,"nan matrix coefficient!");

    // search for "nan" entries in RHS
    //if (prgNAN(RHS,N)) prgError(1,fct,"nan rhs entry!");
  }

  count = 0;
  for(k=0; k<mtx.outerSize(); ++k)
  {
    for(SparseMatrixXd::InnerIterator it(mtx,k); it; ++it)
    {
      //cout << it.row() << '\t' << it.col() << '\t' << it.value() << endl;
      array[count++] = it.value();
    }
  }

  computerTime.go(fct);

  ma41ad_(&JOB,&N,&NE,&row[0], &col[0], &array[0], &rhsVec[0],COLSCA,ROWSCA,KEEP,IS,&MAXIS,S,&MAXS,
               CNTL,ICNTL,INFO,RINFO);
 
  solverTime.total             -= solverTime.factoriseAndSolve;
  solverTime.factoriseAndSolve += computerTime.stop(fct);
  solverTime.total             += solverTime.factoriseAndSolve;

  soln = rhsVec;

  currentStatus = FACTORISE_OK;
  
  if(INFO[0] != 0)
    return 2;
  
  if(checkIO)
  {
    // search for "nan" entries in solution vector
    //if (prgNAN(RHS,N)) prgError(1,fct,"nan entry in solution vector!");
  }

  return 0;
}




void SolverMA41Eigen::free()
{
  if (numSCA > 10) 
  {  
    if (ROWSCA != NULL) delete [] ROWSCA; ROWSCA = new double [10];
    if (COLSCA != NULL) delete [] COLSCA; COLSCA = new double [10];
  }

  if (  IS != NULL) delete [] IS; IS = NULL; MAXIS = 0;

  //mtx.free();

  currentStatus = EMPTY;
  
  return;
}
    



