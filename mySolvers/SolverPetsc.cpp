
#include "Debug.h"
#include "SolverPetsc.h"
#include "petscmat.h"
#include "petscksp.h"
#include "SolverTime.h"
#include "ComputerTime.h"


extern SolverTime      solverTime;
extern ComputerTime    computerTime;


SolverPetsc::SolverPetsc()
{
  if(debug) cout << " SolverPetsc constructor\n\n";

  KSPCreate(PETSC_COMM_WORLD, &ksp);

  //MatCreate(PETSC_COMM_WORLD, &mtx);
  //VecCreate(PETSC_COMM_WORLD, &soln);
  //VecCreate(PETSC_COMM_WORLD, &solnPrev);
  //VecCreate(PETSC_COMM_WORLD, &rhsVec);
  //VecCreate(PETSC_COMM_WORLD, &reac);
}


SolverPetsc::~SolverPetsc()
{
  //if (debug)  cout << " SolverPetsc destructor\n\n";
  //PetscPrintf(MPI_COMM_WORLD, "SolverPetsc::~SolverPetsc() \n");
  //cout << "SolverPetsc::~SolverPetsc() " << endl;
  //free();
  //PetscPrintf(MPI_COMM_WORLD, "SolverPetsc::~SolverPetsc() \n");
  //cout << "SolverPetsc::~SolverPetsc() " << endl;
}


int SolverPetsc::initialise(int p1, int p2, int p3)
{
    nRow = nCol = p3;

    return 0;
}


int SolverPetsc::setSolverAndParameters()
{
    // set the KSP 
    ///////////////////////////////////////////////

    //ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);CHKERRQ(ierr);

    ierr = KSPSetOperators(ksp, mtx, mtx);CHKERRQ(ierr);

    KSPSetType(ksp, KSPCG);

    //KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
    //KSPSetInitialGuessNonzero(ksp, PETSC_FALSE);

    ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

    // set the PC
    ///////////////////////////////////////////////

    //ierr = PCCreate(PETSC_COMM_WORLD, &pc);CHKERRQ(ierr);

    ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);

    ierr = KSPSetReusePreconditioner(ksp, PETSC_FALSE);

    PCSetType(pc, PCILU);
    //PCSetType(pc, PCBJACOBI);

    ierr = PCSetFromOptions(pc);CHKERRQ(ierr);

    //PetscPrintf(MPI_COMM_WORLD, " qqqqqqqqqq \n");

    return 0;
}




int SolverPetsc::zeroMtx()
{
  //cout << " SolverPetsc::zeroMtx() " << endl;
  ierr = MatAssemblyBegin(mtx,MAT_FINAL_ASSEMBLY);//CHKERRQ(ierr);
  ierr = MatAssemblyEnd(mtx,MAT_FINAL_ASSEMBLY);//CHKERRQ(ierr);
  //cout << " SolverPetsc::zeroMtx() " << endl;
  MatZeroEntries(mtx);
  //cout << " SolverPetsc::zeroMtx() " << endl;

  VecAssemblyBegin(rhsVec);
  VecAssemblyEnd(rhsVec);

  VecZeroEntries(rhsVec);

  VecAssemblyBegin(reac);
  VecAssemblyEnd(reac);

  VecZeroEntries(reac);

  return 1;
}



int SolverPetsc::free()
{
  PetscPrintf(MPI_COMM_WORLD, "SolverPetsc::free() \n");

  //ierr = KSPReset(ksp);CHKERRQ(ierr);
  //cout << " ierr = " << ierr << endl;
  ierr = VecDestroy(&soln);CHKERRQ(ierr);
  //cout << " ierr = " << ierr << endl;
  ierr = VecDestroy(&solnPrev);CHKERRQ(ierr);
  //cout << " ierr = " << ierr << endl;
  ierr = VecDestroy(&rhsVec);CHKERRQ(ierr);
  //cout << " ierr = " << ierr << endl;
  ierr = VecDestroy(&reac);CHKERRQ(ierr);
  //cout << " ierr = " << ierr << endl;
  ierr = MatDestroy(&mtx);CHKERRQ(ierr);
  //cout << " ierr = " << ierr << endl;
  //ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  //cout << " ierr = " << ierr << endl;
  //ierr = PCDestroy(&pc);CHKERRQ(ierr);
  //cout << " ierr = " << ierr << endl;
  ierr = PCReset(pc);CHKERRQ(ierr);
  //cout << " ierr = " << ierr << endl;
  //ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  ierr = KSPReset(ksp);CHKERRQ(ierr);
  //cout << " ierr = " << ierr << endl;

  PetscPrintf(MPI_COMM_WORLD, "SolverPetsc::free() \n");
  
  return 1;
}




int SolverPetsc::printInfo()
{
  MatGetInfo(mtx, MAT_LOCAL, &info);

  PetscPrintf(MPI_COMM_WORLD, "Petsc solver:  nRow = %12d \n", nRow);
  PetscPrintf(MPI_COMM_WORLD, "               nnz  = %12d \n\n", info.nz_allocated);

  return 1;
}


int SolverPetsc::printMatrix(int dig, int dig2, bool gfrmt, int indent, bool interactive)
{
  printInfo();

  ierr = MatAssemblyBegin(mtx, MAT_FINAL_ASSEMBLY);//CHKERRQ(ierr);
  ierr = MatAssemblyEnd(mtx, MAT_FINAL_ASSEMBLY);//CHKERRQ(ierr);
 
  //PetscViewerCreate(PETSC_COMM_WORLD, &viewer_matx);
  //PetscViewerDrawOpen();
  //PetscViewerSetFormat(viewer_matx, PETSC_VIEWER_ASCII_MATLAB);
  MatView(mtx, PETSC_VIEWER_STDOUT_WORLD);
  //MatView(A,viewer);

  return 1;
}


double SolverPetsc::giveMatrixCoefficient(int row, int col)
{ 
  //return mtx(row,col,true);
  return 0.0;
}



int SolverPetsc::factorise()
{
  char fct[] = "SolverPetsc::factorise";

  if (currentStatus != ASSEMBLY_OK) { prgWarning(1,fct,"assemble matrix first!"); return 1; }

  if (checkIO)
  {
    // search for "nan" entries in matrix coefficients

    //if (prgNAN(mtx.x.x,NE)) prgError(1,fct,"nan matrix coefficient!");
  }

  computerTime.go(fct);

  currentStatus = FACTORISE_OK;
  
  solverTime.total     -= solverTime.factorise;
  solverTime.factorise += computerTime.stop(fct);
  solverTime.total     += solverTime.factorise;

  return 0;
}


int SolverPetsc::solve()
{
  char fct[] = "SolverPetsc::solve";

  time_t tstart, tend; 

  PetscInt its;
  KSPConvergedReason reason;

  if (currentStatus != FACTORISE_OK) { prgWarning(1,fct,"factorise matrix first!"); return -1; }

  ierr = MatAssemblyBegin(mtx,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(mtx,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  PetscPrintf(MPI_COMM_WORLD, "  SolverPetsc::solve() ... Matrix Assembly ...  \n\n");

  //PetscViewerCreate(PETSC_COMM_WORLD, &viewer_matx);
  //PetscViewerDrawOpen();
  //PetscViewerSetFormat(viewer_matx, PETSC_VIEWER_ASCII_MATLAB);

  //MatView(mtx,PETSC_VIEWER_STDOUT_WORLD);

  ierr = VecAssemblyBegin(rhsVec); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(rhsVec); CHKERRQ(ierr);

  PetscPrintf(MPI_COMM_WORLD, "  SolverPetsc::solve() ... rhsVec Assembly ...  \n\n");

  ierr = VecAssemblyBegin(soln); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(soln); CHKERRQ(ierr);

  PetscPrintf(MPI_COMM_WORLD, "  SolverPetsc::solve() ... soln Assembly ...  \n\n");

  VecZeroEntries(soln);

  PetscPrintf(MPI_COMM_WORLD, "  SolverPetsc::solve() ... vec zero ...  \n\n");

  computerTime.go(fct);

  //VecView(rhsVec, PETSC_VIEWER_STDOUT_WORLD);

  tstart = time(0);

  ierr = KSPSolve(ksp,rhsVec,soln); CHKERRQ(ierr);

  PetscPrintf(MPI_COMM_WORLD, "  SolverPetsc::solve() ... KSP solve ...  \n\n");

  ierr = KSPGetIterationNumber(ksp,&its); CHKERRQ(ierr);

  KSPGetConvergedReason(ksp,&reason);

  if(reason<0)
  {
    PetscPrintf(MPI_COMM_WORLD, "Divergence.\n");
  }
  else
  {
    KSPGetIterationNumber(ksp,&its);
    PetscPrintf(MPI_COMM_WORLD, "Convergence in %d iterations.\n", its);
  }

  //ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);//CHKERRQ(ierr);
  //VecView(soln, PETSC_VIEWER_STDOUT_WORLD);

  tend = time(0); 
  
  //cout << "SolverPetsc::solve()  took "<< difftime(tend, tstart) <<" second(s)."<< endl;

  solverTime.total -= solverTime.solve;
  solverTime.solve += computerTime.stop(fct);
  solverTime.total += solverTime.solve;


  return 0;
}






int SolverPetsc::factoriseAndSolve()
{
  char fct[] = "SolverPetsc::factoriseAndSolve";

  if(currentStatus != ASSEMBLY_OK)
    { prgWarning(1,fct,"assemble matrix first!"); return -1; }

  factorise();
  solve();

  return 0;
}




int SolverPetsc::assembleMatrixAndVector(vector<int>& row, vector<int>& col, MatrixXd& Klocal, VectorXd& Flocal)
{
  int ii, jj;

  int size1 = row.size();
  int size2 = col.size();

  for(ii=0;ii<size1;ii++)
  {
    VecSetValue(rhsVec, row[ii], Flocal(ii), ADD_VALUES);

    for(jj=0;jj<size2;jj++)
    {
      MatSetValue(mtx, row[ii], col[jj], Klocal(ii,jj), ADD_VALUES);
    }
  }

  return 1;
}



int SolverPetsc::assembleMatrixAndVector(int start1, int start2, vector<int>& row, vector<int>& col, MatrixXd& Klocal, VectorXd& Flocal)
{
  int ii, jj, r1, c1;

  int size1 = row.size();
  int size2 = col.size();

  for(ii=0;ii<size1;ii++)
  {
    r1 = start1 + row[ii];

    VecSetValue(rhsVec, r1, Flocal(ii), ADD_VALUES);

    for(jj=0;jj<size2;jj++)
    {
      c1 = start2 + col[jj];
      MatSetValue(mtx, r1, c1, Klocal(ii,jj), ADD_VALUES);
    }
  }

  return 0;
}




int SolverPetsc::assembleMatrixAndVectorCutFEM(int start, int c1, vector<int>& forAssyElem, vector<int>& map_to_cutfem, MatrixXd& Klocal, VectorXd& Flocal)
{
  int ii, jj, r, kk=0;

  int  size1 = forAssyElem.size();
  int  size2 = size1*size1;

  PetscInt  row1[size1];
  PetscScalar  array[size2];

  for(ii=0;ii<size1;ii++)
  {
    r = map_to_cutfem[forAssyElem[ii]];

    VecSetValue(rhsVec, r, Flocal(ii), ADD_VALUES);

    row1[ii] = r;

    for(jj=0;jj<size1;jj++)
    {
      array[kk++] = Klocal(ii,jj);
    }
  }

  MatSetValues(mtx, size1, row1, size1, row1, array, ADD_VALUES);

  return 0;
}



int SolverPetsc::assembleMatrixAndVectorCutFEM2(int start1, int start2, vector<int>& tempVec, 
     vector<int>& forAssy1, vector<int>& forAssy2, MatrixXd& Klocal, VectorXd& Flocal1, VectorXd& Flocal2)
{
  // to assemble coupling matrices arriving from 
  // jump conditions in the gradients across an interface

  int ii, jj, size1, kk, r, c;

  //printVector(vec1);
  //printVector(vec2);
  
  size1 = tempVec.size();

  for(ii=0;ii<size1;ii++)
  {
    // force vector for domain #2
    r = start2 + forAssy2[tempVec[ii]];

    VecSetValues(rhsVec, 1, &r, &Flocal2(ii), ADD_VALUES);

    // force vector for domain #1
    r = start1 + forAssy1[tempVec[ii]];

    VecSetValues(rhsVec, 1, &r, &Flocal1(ii), ADD_VALUES);

    // coupling matrix
    for(jj=0;jj<size1;jj++)
    {
      c = start2 + forAssy2[tempVec[jj]];
      MatSetValues(mtx, 1, &r, 1, &c, &(Klocal(ii,jj)), ADD_VALUES);
      MatSetValues(mtx, 1, &c, 1, &r, &(Klocal(ii,jj)), ADD_VALUES);
    }
  }

  return 0;
}



int SolverPetsc::assembleMatrixAndVectorCutFEM3(int start1, int start2, vector<int>& row, vector<int>& col, 
     vector<int>& forAssy1, vector<int>& forAssy2, MatrixXd& Klocal, VectorXd& Flocal1)
{
 
  return 0;
}





int SolverPetsc::assembleMatrixAndVectorMixedFormulation(int start, int c1, vector<int>& vec1, vector<int>& vec2, MatrixXd& Klocal, VectorXd& Flocal)
{
  int ii, jj, aa, bb, size1, size2;

  //printVector(vec1);
  //printVector(vec2);
  
  size1 = vec1.size();
  size2 = vec2.size();

  for(ii=0;ii<size1;ii++)
  {
    VecSetValues(rhsVec, 1, &vec1[ii], &Flocal(ii), ADD_VALUES);

    //for(jj=0;jj<size1;jj++)
      //mtx.coeffRef(vec1[ii], vec1[jj]) += Klocal(ii, jj);

    for(jj=0;jj<size2;jj++)
    {
      aa = start + vec2[jj];
      bb = size1 + jj;
      
      MatSetValues(mtx, 1, &vec1[ii], 1, &aa,       &(Klocal(ii,bb)), ADD_VALUES);
      MatSetValues(mtx, 1, &aa,       1, &vec1[ii], &(Klocal(bb,ii)), ADD_VALUES);
    }
  }

  for(ii=0;ii<size2;ii++)
  {
    aa = start + vec2[ii];
    VecSetValues(rhsVec, 1, &aa, &Flocal(size1+ii), ADD_VALUES);
  }
  return 0;
}



int SolverPetsc::assembleVector(int start, int c1, vector<int>& vec1, VectorXd& Flocal)
{
  int ii, jj;

  for(ii=0;ii<vec1.size();ii++)
  {
    VecSetValues(rhsVec, 1, &vec1[ii], &Flocal(ii), ADD_VALUES);
  }

  return 0;
}





