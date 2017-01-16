
#include <iostream>

#include "SolverMA41.h"
#include "FunctionsSolver.h"
#include "FunctionsProgram.h"
#include "SolverWorkSpace.h"
#include "SolverTime.h"
#include "ComputerTime.h"


extern SolverWorkSpace solverWorkSpace;
extern SolverTime      solverTime;
extern ComputerTime    computerTime;


using namespace std;




SolverMA41::SolverMA41(void)
{
  //comprMtxFlg = true;

  ROWSCA = new double [10];
  COLSCA = new double [10];

  IS     = NULL;
  
  numSCA = 10;

  return;
}




SolverMA41::~SolverMA41()
{
  free();

  if (ROWSCA != NULL) delete [] ROWSCA;
  if (COLSCA != NULL) delete [] COLSCA;
  if (IS     != NULL) delete [] IS;

  return;
}





int SolverMA41::initialise(int, int, int)
{
  if (currentStatus != PATTERN_OK)
    { prgWarning(1,"SolverMA41::initialise","prepare matrix pattern first!"); return 1; }

  int    *tmp, i, j, m, 
	 MAXS, JOB, 
	 N  = mtx.nRow, 
	 NE = mtx.x.n, 
	 *IRN, *JCN;

  double *S, RHS[10], *ASPK;
  
  if (mtx.nRow != mtx.nCol) prgError(1,"SolverMA41::initialise","nRow != nCol");

  // initialise controls
	
  ma41id_(CNTL,ICNTL,KEEP);

  CNTL [0] = 0.01;//.01;  // pivoting
//  ICNTL[0] = 1;
  //ICNTL[7] = 6;    // scaling
  //ICNTL[9] = 10;   // max iterations


//    for(int ii=0;ii<10;ii++)       cout << " ICNTL[" << ii << "] .... : " << ICNTL[ii] << endl;  cout << endl;

//  ICNTL[4] = 3;
  
  // analyse matrix pattern

  MAXIS = 2 * NE + 11*N + 1;
  
  IS    = new int [MAXIS];
  
  JOB   = 1;

  MAXS  = solverWorkSpace.dimDbl;

  S     = solverWorkSpace.dbl;
  
  IRN   = mtx.row.x;
  JCN   = mtx.col.x;

  ASPK  = mtx.x.x;
  
  ma41ad_(&JOB,&N,&NE,IRN,JCN,ASPK,RHS,COLSCA,ROWSCA,KEEP,IS,&MAXIS,S,&MAXS,
		  CNTL,ICNTL,INFO,RINFO);

  // allocate memory for IS if required

//  m = INFO[6] + 10 + roundToInt(0.2 * (double)NE);

  m = roundToInt(2.0 * (double)INFO[6]);
  
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
  
  m = roundToInt(3.0 * (double)INFO[7]);
  
//  cout << " m .... : " << m << endl;
//  cout << " KEEP[12] .... : " << KEEP[12] << endl;
  
  solverWorkSpace.expand(m);

  currentStatus = INIT_OK;

  return 0;
}





int SolverMA41::factorise(void)
{
  char fct[] = "SolverMA41::factorise";

  if (currentStatus != ASSEMBLY_OK) { prgWarning(1,fct,"assemble matrix first!"); return 1; }

  int    JOB  = 2, 
         MAXS = solverWorkSpace.dimDbl,
         *IRN = mtx.row.x,
         *JCN = mtx.col.x,
         N    = mtx.nRow,
         NE   = mtx.x.n;

  double *ASPK = mtx.x.x,
	 *S = solverWorkSpace.dbl,
	 RHS[10];

  if (checkIO)
  {
    // search for "nan" entries in matrix coefficients

    if (prgNAN(mtx.x.x,NE)) prgError(1,fct,"nan matrix coefficient!");
  }

  computerTime.go(fct);

  ma41ad_(&JOB,&N,&NE,IRN,JCN,ASPK,RHS,COLSCA,ROWSCA,KEEP,IS,&MAXIS,S,&MAXS,
               CNTL,ICNTL,INFO,RINFO);
 
  currentStatus = FACTORISE_OK;
  
  solverTime.total     -= solverTime.factorise;
  solverTime.factorise += computerTime.stop(fct);
  solverTime.total     += solverTime.factorise;

  return INFO[0];
}






double *SolverMA41::solve(double *RHS, int)
{
  char fct[] = "SolverMA41::solve";

  if (currentStatus != FACTORISE_OK) { prgWarning(1,fct,"factorise matrix first!"); return NULL; }

  int    JOB  = 3, 
         MAXS = solverWorkSpace.dimDbl,
         *IRN = mtx.row.x,
         *JCN = mtx.col.x,
         N    = mtx.nRow,
         NE   = mtx.x.n;

  double *ASPK = mtx.x.x,
	 *S = solverWorkSpace.dbl;
  
  if (checkIO)
  {
    // search for "nan" entries in RHS

    if (prgNAN(RHS,N)) prgError(1,fct,"nan rhs entry!");
  }

  computerTime.go(fct);

  ma41ad_(&JOB,&N,&NE,IRN,JCN,ASPK,RHS,COLSCA,ROWSCA,KEEP,IS,&MAXIS,S,&MAXS,
             CNTL,ICNTL,INFO,RINFO);

  solverTime.total -= solverTime.solve;
  solverTime.solve += computerTime.stop(fct);
  solverTime.total += solverTime.solve;

  if (INFO[0] != 0) return NULL;
  
  if (checkIO)
  {
    // search for "nan" entries in solution vector

    if (prgNAN(RHS,N)) prgError(1,fct,"nan entry in solution vector!");
  }

  return RHS;
}






double *SolverMA41::factoriseAndSolve(double *RHS, int)
{
  char fct[] = "SolverMA41::factoriseAndSolve";

  if (currentStatus != ASSEMBLY_OK) { prgWarning(1,fct,"assemble matrix first!"); return NULL; }

  int    JOB  = 5, 
         MAXS = solverWorkSpace.dimDbl,
         *IRN = mtx.row.x,
         *JCN = mtx.col.x,
         N    = mtx.nRow,
         NE   = mtx.x.n;

  double *ASPK = mtx.x.x,
	 *S = solverWorkSpace.dbl;

  if (checkIO)
  {
    // search for "nan" entries in matrix coefficients

    if (prgNAN(mtx.x.x,NE)) prgError(1,fct,"nan matrix coefficient!");

    // search for "nan" entries in RHS

    if (prgNAN(RHS,N)) prgError(1,fct,"nan rhs entry!");
  }


  computerTime.go(fct);

  ma41ad_(&JOB,&N,&NE,IRN,JCN,ASPK,RHS,COLSCA,ROWSCA,KEEP,IS,&MAXIS,S,&MAXS,
               CNTL,ICNTL,INFO,RINFO);
 
  solverTime.total             -= solverTime.factoriseAndSolve;
  solverTime.factoriseAndSolve += computerTime.stop(fct);
  solverTime.total             += solverTime.factoriseAndSolve;

  currentStatus = FACTORISE_OK;
  
  if (INFO[0] != 0) return NULL;
  
  if (checkIO)
  {
    // search for "nan" entries in solution vector

    if (prgNAN(RHS,N)) prgError(1,fct,"nan entry in solution vector!");
  }

  return RHS;
}







void SolverMA41::free(void)
{
  if (numSCA > 10) 
  {  
    if (ROWSCA != NULL) delete [] ROWSCA; ROWSCA = new double [10];
    if (COLSCA != NULL) delete [] COLSCA; COLSCA = new double [10];
  }

  if (  IS != NULL) delete [] IS; IS = NULL; MAXIS = 0;

  mtx.free();

  currentStatus = EMPTY;
  
  return;
}
    



