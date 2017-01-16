
#include <iostream>

#include "SolverPardisoPetsc.h"
#include "FunctionsSolver.h"
#include "FunctionsProgram.h"
#include "SolverTime.h"
#include "ComputerTime.h"
#include "util.h"

#include "petscksp.h"
#include "petscmat.h"
#include "petscsys.h"


extern SolverTime      solverTime;
extern ComputerTime    computerTime;
extern int             countThis;

using namespace std;




SolverPardisoPetsc::SolverPardisoPetsc()
{
  comprMtxFlg = true;

  return;
}




SolverPardisoPetsc::~SolverPardisoPetsc()
{
  int phase = -1, error = 0, idum;

  double ddum;

  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &N, array, csr, col, &idum, &NRHS,
           IPARM, &MSGLVL, &ddum, &ddum, &error, DPARM);

  return;
}





int SolverPardisoPetsc::initialise(int numProc, int matrixType, int nrow)
{
  cout << " ppppppppppp " << endl;

  SolverPetsc::initialise(numProc, matrixType, nrow);

  cout << " ppppppppppp " << endl;

  char fct[] = "SolverPardisoPetsc::initialise";

  //cout << fct << ": " << ++countThis << "\n\n";

  //if (currentStatus != PATTERN_OK)
    //{ prgWarning(1,fct,"prepare matrix pattern first!"); return 1; }

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

  N = nrow; // dimension of matrix

  nx = N;

  //if(X != NULL)
    //delete [] X;
  //X = NULL;
 
  X = new double [nx];

  pardisoinit_(PT, &MTYPE, &SOLVER, IPARM, DPARM, &error);

  if (error != 0)
  {
    if (error == -10) prgError(1,fct,"no license file found.");
    if (error == -11) prgError(2,fct,"license is expired.");
    if (error == -12) prgError(3,fct,"wrong username or hostname.");
  }

  COUT << "\n\n PARDISO license check was successful.\n\n";

  ierr = MatAssemblyBegin(mtx,MAT_FINAL_ASSEMBLY);//CHKERRQ(ierr);
  ierr = MatAssemblyEnd(mtx,MAT_FINAL_ASSEMBLY);//CHKERRQ(ierr);
  
  //MatView(mtx, PETSC_VIEWER_STDOUT_WORLD);

  //cout << " ppppppppppp " << endl;

  PetscBool flag;

  //ierr = MatGetRowIJ(mtx, 1, PETSC_FALSE, PETSC_FALSE, &nRow, &csr, &col, &flag);

  //cout << " llllllllll " << endl;

  ierr = MatSeqAIJGetArray(mtx, &array);
  //ierr = MatGetArray(mtx, &array);

  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &N, array, csr, col, &idum, &NRHS,
           IPARM, &MSGLVL, &ddum, &ddum, &error, DPARM);

  cout << " llllllllll " << endl;

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







int SolverPardisoPetsc::factorise(void)
{
  char fct[] = "SolverPARDISO::factorise";

  if (currentStatus != ASSEMBLY_OK) 
    { prgWarning(1,fct,"assemble matrix first!"); return 1; }

  int phase = 22, idum, error = 0;

  double ddum;

  computerTime.go(fct);

  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &N, array, csr, col, &idum, &NRHS,
           IPARM, &MSGLVL, &ddum, &ddum, &error, DPARM);

  solverTime.total     -= solverTime.factorise;
  solverTime.factorise += computerTime.stop(fct);
  solverTime.total     += solverTime.factorise;

  currentStatus = FACTORISE_OK;
  
  return 0;
}






double *SolverPardisoPetsc::solve(double *RHS, int nrhs)
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
           &N, array, csr, col, &idum, &NRHS,
           IPARM, &MSGLVL, RHS, X, &error, DPARM);
 
  solverTime.total -= solverTime.solve;
  solverTime.solve += computerTime.stop(fct);
  solverTime.total += solverTime.solve;

  return RHS;
}






double *SolverPardisoPetsc::factoriseAndSolve(double *RHS, int nrhs)
{
  char fct[] = "SolverPARDISO::factoriseAndSolve";

  if (currentStatus != ASSEMBLY_OK) 
    { prgWarning(1,fct,"assemble matrix first!"); return NULL; }

  int phase = 23, idum, error = 0;

  N = nRow; // dimension of matrix

  nx = N;

  NRHS = nrhs;
  
  //cout << nx << '\t' << NRHS << endl;

  if (nx < NRHS * N) 
  {
    delete [] X; 
    nx = NRHS * N;
    X = new double [nx];
  }

  ierr = MatAssemblyBegin(mtx,MAT_FINAL_ASSEMBLY);//CHKERRQ(ierr);
  ierr = MatAssemblyEnd(mtx,MAT_FINAL_ASSEMBLY);//CHKERRQ(ierr);
  
  //MatView(mtx, PETSC_VIEWER_STDOUT_WORLD);

  //cout << fct << " NRHS = " << NRHS << endl;

  //cout << " rhsVec - pardiso " << endl;        printVector(RHS, nRow);

  computerTime.go(fct);

  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &N, array, csr, col, &idum, &NRHS,
           IPARM, &MSGLVL, RHS, X, &error, DPARM);

  solverTime.total             -= solverTime.factoriseAndSolve;
  solverTime.factoriseAndSolve += computerTime.stop(fct);
  solverTime.total             += solverTime.factoriseAndSolve;

  currentStatus = FACTORISE_OK;

  //cout << " soln - pardiso " << endl;        printVector(RHS, nRow);

  //cout << " ppppppppppp " << endl;
  
  return RHS;
}







int SolverPardisoPetsc::free()
{
  int phase = -1, idum, error = 0;

  double ddum;

  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &N, array, csr, col, &idum, &NRHS,
           IPARM, &MSGLVL, &ddum, &ddum, &error, DPARM);


  delete [] X;
  X = NULL;

  //mtx.free();
  //compr.free();

  currentStatus = EMPTY;

  return 1;
}
    



