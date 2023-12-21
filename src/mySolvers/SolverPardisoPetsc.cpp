
#include <iostream>

#include "SolverPardisoPetsc.h"
#include "pardiso.h"
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
  //MatCreate(PETSC_COMM_WORLD, &mtx);

  //VecCreate(PETSC_COMM_WORLD, &soln);
  //VecCreate(PETSC_COMM_WORLD, &solnPrev);
  //VecCreate(PETSC_COMM_WORLD, &rhsVec);
  //VecCreate(PETSC_COMM_WORLD, &reac);

#ifndef PARDISO_SOLVER
  cerr << "PARDISO solver is not available ... " << endl;
  cerr << " Need to be compiled with PARDISO_SOLVER preprocessor directive " << endl;
  exit(1);
#endif
}




SolverPardisoPetsc::~SolverPardisoPetsc()
{
  //free();

  phase = -1; error = 0;

#ifdef PARDISO_SOLVER
  pardiso(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &nRow, array, csr, col, perm, &NRHS,
           IPARM, &MSGLVL, &ddum, &ddum, &error, DPARM);
#endif
}




int SolverPardisoPetsc::initialise(int numProc, int matrixType, int nrow1)
{
  nRow = nCol = nrow1;

  rhsTemp.resize(nRow, 0.0);
  solnTemp.resize(nRow, 0.0);

  char fct[] = "SolverPardisoPetsc::initialise";

  if (currentStatus != PATTERN_OK)
  {
    prgWarning(1,fct,"prepare matrix pattern first!");
    return 2;
  }

  phase = 11; error = 0;

  int idum, mtxType = 11;

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

  tmp = getenv("OMP_NUM_THREADS");

  if (tmp != NULL) 
  {
    sscanf(tmp,"%d", &idum);
    //if (idum != IPARM[2]) prgError(1,fct,"set environment variable OMP_NUM_THREADS to numProc!");
  }
  else prgError(2,fct,"set environment variable OMP_NUM_THREADS!");

  IPARM[2] = max(1,idum);  // number of processors (no default value available)


#ifdef PARDISO_SOLVER
  pardisoinit(PT, &MTYPE, &SOLVER, IPARM, DPARM, &error);
#endif

  if (error != 0)
  {
    if (error == -10) prgError(1,fct,"no license file found.");
    if (error == -11) prgError(2,fct,"license is expired.");
    if (error == -12) prgError(3,fct,"wrong username or hostname.");
  }

  COUT << "\n\n PARDISO license check was successful.\n\n";

  ierr = MatAssemblyBegin(mtx,MAT_FINAL_ASSEMBLY);//CHKERRQ(ierr);
  ierr = MatAssemblyEnd(mtx,MAT_FINAL_ASSEMBLY);//CHKERRQ(ierr);

  PetscBool flag;

  ierr = MatGetRowIJ(mtx, 1, PETSC_FALSE, PETSC_FALSE, &nRow, &csr, &col, &flag);

  ierr = MatSeqAIJGetArray(mtx, &array);
  //ierr = MatGetArray(mtx, &array);

#ifdef PARDISO_SOLVER
  pardiso(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &nRow, array, csr, col, &idum, &NRHS,
           IPARM, &MSGLVL, &ddum, &ddum, &error, DPARM);
#endif
  if (error != 0)
  {
    COUT << "PARDISO ERROR = " << error << "\n\n";
    prgError(4,fct,"symbolic factorisation failed.");
  }

  //IPARM[5] = 1; // overwrite RHS with solution
  IPARM[5] = 0; // do not overwrite RHS with solution
  IPARM[7] = 1; // max number of iterative refinement steps

  currentStatus = INIT_OK;

  return 0;
}



int SolverPardisoPetsc::factorise()
{
  char fct[] = "SolverPARDISO::factorise";

  if (currentStatus != ASSEMBLY_OK) 
    { prgWarning(1,fct,"assemble matrix first!"); return 1; }

  phase = 22; error = 0;

  computerTime.go(fct);

#ifdef PARDISO_SOLVER
  pardiso(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &nRow, array, csr, col, perm, &NRHS,
           IPARM, &MSGLVL, &ddum, &ddum, &error, DPARM);
#endif

  solverTime.total     -= solverTime.factorise;
  solverTime.factorise += computerTime.stop(fct);
  solverTime.total     += solverTime.factorise;

  currentStatus = FACTORISE_OK;

  return 0;
}






int  SolverPardisoPetsc::solve()
{ 
  char fct[] = "SolverPARDISO::solve";
 
  if(currentStatus != FACTORISE_OK)
    { prgWarning(1,fct,"factorise matrix first!"); return 0; }

  phase = 33; error = 0;

  computerTime.go(fct);

  PetscScalar *arrayRhs, *arraySoln;

  VecGetArray(rhsVec, &arrayRhs);
  VecGetArray(soln, &arraySoln);

#ifdef PARDISO_SOLVER
  pardiso(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &nRow, array, csr, col, perm, &NRHS,
           IPARM, &MSGLVL, arrayRhs, arraySoln, &error, DPARM);
#endif

  //for(int ii=0; ii<nRow; ii++)
    //VecSetValue(soln, ii, arraySoln[ii], INSERT_VALUES);

  //VecAssemblyBegin(soln);
  //VecAssemblyEnd(soln);

  VecRestoreArray(rhsVec, &arrayRhs);
  VecRestoreArray(soln, &arraySoln);

  solverTime.total -= solverTime.solve;
  solverTime.solve += computerTime.stop(fct);
  solverTime.total += solverTime.solve;

  return 0;
}






int  SolverPardisoPetsc::factoriseAndSolve()
{
  char fct[] = "SolverPARDISO::factoriseAndSolve";

  if (currentStatus != ASSEMBLY_OK) 
    { prgWarning(1,fct,"assemble matrix first!"); return 0; }

  phase = 23; error = 0;

  ierr = MatAssemblyBegin(mtx,MAT_FINAL_ASSEMBLY);//CHKERRQ(ierr);
  ierr = MatAssemblyEnd(mtx,MAT_FINAL_ASSEMBLY);//CHKERRQ(ierr);

  //MatView(mtx, PETSC_VIEWER_STDOUT_WORLD);

  VecAssemblyBegin(rhsVec);
  VecAssemblyEnd(rhsVec);

  computerTime.go(fct);

  PetscScalar *arrayRhs;
  VecGetArray(rhsVec, &arrayRhs);

  int ii=0;
  for(ii=0; ii<nRow; ii++)
  {
    rhsTemp[ii] = arrayRhs[ii];
    solnTemp[ii] = 0.0;
  }

  VecRestoreArray(rhsVec, &arrayRhs);

#ifdef PARDISO_SOLVER
  pardiso(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &nRow, array, csr, col, perm, &NRHS,
           IPARM, &MSGLVL, &rhsTemp[0], &solnTemp[0], &error, DPARM);
#endif

  VecZeroEntries(soln);
  for(ii=0; ii<nRow; ii++)
  {
    //cout << ii << '\t' << rhsTemp[ii] << '\t' << solnTemp[ii] << endl;
    ierr = VecSetValue(soln, ii, solnTemp[ii], ADD_VALUES);
  }

  VecAssemblyBegin(soln);
  VecAssemblyEnd(soln);

  //printf("Peak memory [kB] phase 1       = %d \n", IPARM[14]);
  //printf("Permanent integer memory [kb]. = %d \n", IPARM[15]);
  //printf("Peak real memory [kB]          = %d \n", IPARM[16]);
  //printf("Number of nonzeros in LU.      = %d \n", IPARM[17]);
  //printf("Gflops for LU factorization.   = %d \n", IPARM[18]);

  solverTime.total             -= solverTime.factoriseAndSolve;
  solverTime.factoriseAndSolve += computerTime.stop(fct);
  solverTime.total             += solverTime.factoriseAndSolve;

  currentStatus = FACTORISE_OK;

  return 0;
}





int SolverPardisoPetsc::free()
{
  ierr = VecDestroy(&soln);CHKERRQ(ierr);
  ierr = VecDestroy(&solnPrev);CHKERRQ(ierr);
  ierr = VecDestroy(&rhsVec);CHKERRQ(ierr);
  ierr = VecDestroy(&reac);CHKERRQ(ierr);

  phase = -1; error = 0;

#ifdef PARDISO_SOLVER
  pardiso(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &nRow, array, csr, col, perm, &NRHS,
           IPARM, &MSGLVL, &ddum, &ddum, &error, DPARM);
#endif

  PetscBool flag;

  ierr = MatRestoreRowIJ(mtx, 1, PETSC_FALSE, PETSC_FALSE, &nRow, &csr, &col, &flag);

  MatSeqAIJRestoreArray(mtx, &array);

  ierr = MatDestroy(&mtx);CHKERRQ(ierr);
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);

  currentStatus = EMPTY;

  return 0;
}



