
#include <iostream>

#include "SolverPARDISO.h"
#include "FunctionsSolver.h"
#include "FunctionsProgram.h"
#include "SolverTime.h"
#include "ComputerTime.h"


extern SolverWorkSpace solverWorkSpace;
extern SolverTime      solverTime;
extern ComputerTime    computerTime;
extern int             countThis;

using namespace std;




SolverPARDISO::SolverPARDISO(void)
{
  comprMtxFlg = true;

  return;
}




SolverPARDISO::~SolverPARDISO()
{
  int phase = -1, error = 0, idum;

  double ddum;

  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &N, mtx.x.x, compr.x, mtx.col.x, &idum, &NRHS,
           IPARM, &MSGLVL, &ddum, &ddum, &error, DPARM);

  return;
}





int SolverPARDISO::initialise(int numProc, int matrixType, int)
{
  char fct[] = "SolverPARDISO::initialise";

  //cout << fct << ": " << ++countThis << "\n\n";

  if (currentStatus != PATTERN_OK)
    { prgWarning(1,fct,"prepare matrix pattern first!"); return 1; }

  double ddum;

  int idum, phase = 11, error = 0, mtxType = 11;

  char *tmp;

  switch (matrixType)
  {
    case PARDISO_STRUCT_SYM: mtxType =  1; break; // real and structurally symmetric

    case PARDISO_UNSYM:      mtxType = 11; break; // real and unsymmetric

    default:                 prgWarning(1,fct,"matrix assumed to be unsymmetric!");
  }

  SOLVER = 0;       // sparse direct solver
  MTYPE  = mtxType; // matrix type
  MAXFCT = 1;       // maximum number of factorisations of same sparsity pattern
                    //      to be kept in memory
  MNUM   = 1;       // which factorisation to use
  NRHS   = 1;       // number of right hand sides
  MSGLVL = 0;       // output message level (1 -> print statistical information)

  IPARM[0] = 0;     // PARADISO will set IPARM to default values

  IPARM[2] = max(1,numProc);  // number of processors (no default value available)
  
  //cout << " numProc " << numProc << endl;
  
  //IPARM[9] = 1;

  tmp = getenv("OMP_NUM_THREADS");

  if (tmp != NULL) 
  {
    sscanf(tmp,"%d", &idum);
    if (idum != IPARM[2]) prgError(1,fct,"set environment variable OMP_NUM_THREADS to numProc!");
  }
  else prgError(2,fct,"set environment variable OMP_NUM_THREADS!");

  N = mtx.nRow, // dimension of matrix

  nx = N;

  //if(X != NULL)
    //delete [] X;
  //X = NULL;
 
  X = new double [nx];

  pardisoinit_(PT,&MTYPE,&SOLVER,IPARM,DPARM,&error);

  if (error != 0)
  {
    if (error == -10) prgError(1,fct,"no license file found.");
    if (error == -11) prgError(2,fct,"license is expired.");
    if (error == -12) prgError(3,fct,"wrong username or hostname.");
  }
  //COUT << "PARDISO license check was successful.\n\n";

  //cout << N << "," << compr[N]-1 << "\n";

  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &N, mtx.x.x, compr.x, mtx.col.x, &idum, &NRHS,
           IPARM, &MSGLVL, &ddum, &ddum, &error, DPARM);

  if (error != 0)
  {
    COUT << "PARDISO ERROR = " << error << "\n\n";
    prgError(4,fct,"symbolic factorisation failed.");
  }

  IPARM[5] = 1; // overwrite RHS with solution
  IPARM[7] = 1; // max number of iterative refinement steps

  currentStatus = INIT_OK;

  return 0;
}







int SolverPARDISO::factorise(void)
{
  char fct[] = "SolverPARDISO::factorise";

  if (currentStatus != ASSEMBLY_OK) 
    { prgWarning(1,fct,"assemble matrix first!"); return 1; }

  int phase = 22, idum, error = 0;

  double ddum;

  computerTime.go(fct);

  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &N, mtx.x.x, compr.x, mtx.col.x, &idum, &NRHS,
           IPARM, &MSGLVL, &ddum, &ddum, &error, DPARM);

  solverTime.total     -= solverTime.factorise;
  solverTime.factorise += computerTime.stop(fct);
  solverTime.total     += solverTime.factorise;

  currentStatus = FACTORISE_OK;
  
  return 0;
}






double *SolverPARDISO::solve(double *RHS, int nrhs)
{ 
  char fct[] = "SolverPARDISO::solve";
 
  if (currentStatus != FACTORISE_OK)
    { prgWarning(1,fct,"factorise matrix first!"); return NULL; }

  int phase = 33, idum, error = 0;

  NRHS = nrhs;
  
  //cout << nx << '\t' << NRHS << endl;

  if (nx < NRHS * N) 
  {
    delete [] X; 
    nx = NRHS * N;
    X = new double [nx];
  }

  //cout << fct << " NRHS = " << NRHS << "\n";

  computerTime.go(fct);

  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &N, mtx.x.x, compr.x, mtx.col.x, &idum, &NRHS,
           IPARM, &MSGLVL, RHS, X, &error, DPARM);
 
  solverTime.total -= solverTime.solve;
  solverTime.solve += computerTime.stop(fct);
  solverTime.total += solverTime.solve;

  return RHS;
}






double *SolverPARDISO::factoriseAndSolve(double *RHS, int nrhs)
{
  if(mtx.nRow < 20)
  {
    cout << " CSRrow ... " << endl;
    cout << compr << endl;
    cout << endl;

    cout << " CSRcol ... " << endl;
    cout << mtx.col << endl;
    cout << endl;

    cout << " CSRval ... " << endl;
    cout << mtx.x << endl;
    cout << endl;
  }

  char fct[] = "SolverPARDISO::factoriseAndSolve";

  if (currentStatus != ASSEMBLY_OK) 
    { prgWarning(1,fct,"assemble matrix first!"); return NULL; }

  int phase = 23, idum, error = 0;

  NRHS = nrhs;
  
  cout << nx << '\t' << NRHS << endl;

  if (nx < NRHS * N) 
  {
    delete [] X; 
    nx = NRHS * N;
    X = new double [nx];
  }

  //cout << fct << " NRHS = " << NRHS << "\n";

  computerTime.go(fct);

  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &N, mtx.x.x, compr.x, mtx.col.x, &idum, &NRHS,
           IPARM, &MSGLVL, RHS, X, &error, DPARM);

  solverTime.total             -= solverTime.factoriseAndSolve;
  solverTime.factoriseAndSolve += computerTime.stop(fct);
  solverTime.total             += solverTime.factoriseAndSolve;

  currentStatus = FACTORISE_OK;
  
  return RHS;
}







void SolverPARDISO::free(void)
{
  int phase = -1, idum, error = 0;

  double ddum;

  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &N, mtx.x.x, compr.x, mtx.col.x, &idum, &NRHS,
           IPARM, &MSGLVL, &ddum, &ddum, &error, DPARM);


  delete [] X;
  X = NULL;

  mtx.free();
  compr.free();

  currentStatus = EMPTY;

  return;
}
    



