
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

  //PetscInitialize(NULL,NULL,(char *)0,NULL);

  //PetscInitialize(NULL, NULL, "petsc_options.dat", NULL);

  //MPI_Comm_size(PETSC_COMM_WORLD, &size);

  //if (size != 1) SETERRQ(1,"This is a uniprocessor example only!");

  //ierr = PetscOptionsGetInt(PETSC_NULL, "-n", &nRow, PETSC_NULL);CHKERRQ(ierr);

  //MatCreate(PETSC_COMM_WORLD, &mtx);

}


SolverPetsc::~SolverPetsc()
{
  if (debug)  cout << " SolverPetsc destructor\n\n";

  free();

  //ierr = PetscFinalize();//CHKERRQ(ierr);

  //return;
}


int SolverPetsc::initialise(int p1, int p2, int p3)
{
    nRow = nCol = p3;

    //cout << " p1   = " << p1 << endl;
    //cout << " nRow = " << nRow << endl;

    //VecCreate(PETSC_COMM_WORLD, &soln);
    //VecCreate(PETSC_COMM_WORLD, &solnPrev);
    //VecCreate(PETSC_COMM_WORLD, &rhsVec);
    //VecCreate(PETSC_COMM_WORLD, &reac);
    ////ierr = PetscObjectSetName((PetscObject) x, "Solution");CHKERRQ(ierr);

    //ierr = VecSetSizes(soln, PETSC_DECIDE, nRow);CHKERRQ(ierr);
    //ierr = VecSetSizes(solnPrev, PETSC_DECIDE, nRow);CHKERRQ(ierr);
    //ierr = VecSetSizes(rhsVec, PETSC_DECIDE, nRow);CHKERRQ(ierr);
    //ierr = VecSetSizes(reac, PETSC_DECIDE, p1);CHKERRQ(ierr);

    //ierr = VecSetFromOptions(soln);CHKERRQ(ierr);
    //ierr = VecDuplicate(soln, &rhsVec);CHKERRQ(ierr);
    //ierr = VecDuplicate(soln, &solnPrev);CHKERRQ(ierr);
    //ierr = VecSetFromOptions(reac);CHKERRQ(ierr);

    //cout << " AAAAAAAAAA " << endl;
    //ierr = MatCreateSeqAIJ(PETSC_COMM_WORLD,nRow,nRow,15000,NULL,&mtx);CHKERRQ(ierr);

    //ierr = MatSetSizes(mtx,PETSC_DECIDE,PETSC_DECIDE,nRow,nRow);CHKERRQ(ierr);

    //ierr = MatSetFromOptions(mtx);CHKERRQ(ierr);
    
    //cout << " AAAAAAAAAA " << endl;
    //SetSolverAndParameters();

    return 0;
}


int SolverPetsc::SetSolverAndParameters()
{
    ///////////////////////
    // Create the linear solver and set various options
    ///////////////////////
    PetscPrintf(MPI_COMM_WORLD, " AAAAAAAAAA \n");
    ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);CHKERRQ(ierr);

    ///////////////////////
    //Set operators. Here the matrix that defines the linear system
    //also serves as the preconditioning matrix.
    ///////////////////////
    PetscPrintf(MPI_COMM_WORLD, " BBBBBBBBBB \n");
    ierr = KSPSetOperators(ksp, mtx, mtx);CHKERRQ(ierr);
    PetscPrintf(MPI_COMM_WORLD, " AAAAAAAAAA \n");
    ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
    PetscPrintf(MPI_COMM_WORLD, " cccccccccc \n");
    ///////////////////////
    /*
    Set linear solver defaults for this problem (optional).
    - By extracting the KSP and PC contexts from the KSP context,
    we can then directly call any KSP and PC routines to set
    various options.
    - The following four statements are optional; all of these
    parameters could alternatively be specified at runtime via
    KSPSetFromOptions();
    */
    ///////////////////////

    //KSPSetType(ksp, KSPPREONLY);
    KSPSetType(ksp, KSPCG);
    //KSPSetType(ksp, KSPGMRES);
    //KSPSetType(ksp, KSPBCGS);
    //KSPSetType(ksp, KSPBICG);
    //KSPSetType(ksp, KSPLSQR);
    //KSPGMRESSetRestart(ksp, 100);
    PetscPrintf(MPI_COMM_WORLD, " AAAAAAAAAA \n");

    ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);

    //ierr = PCSetType(pc,PCJACOBI);CHKERRQ(ierr);
    ierr = PCSetType(pc,PCILU);CHKERRQ(ierr);
    //ierr = PCSetType(pc,PCLU);CHKERRQ(ierr);
    //ierr = PCSetType(pc,PCCholesky);CHKERRQ(ierr);
    //ierr = PCSetType(pc,PCNONE);CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp,1.0e-12,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
    PetscPrintf(MPI_COMM_WORLD, " qqqqqqqqqq \n");

    // Set runtime options, e.g.,
    // -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
    // These options will override those specified above as long as
    // KSPSetFromOptions() is called _after_ any other customization
    // routines.

    //KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
    KSPSetInitialGuessNonzero(ksp, PETSC_FALSE);

    ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

    PetscPrintf(MPI_COMM_WORLD, " qqqqqqqqqq \n");

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
    //ierr = VecDestroy(&soln);CHKERRQ(ierr);
    //ierr = VecDestroy(&solnPrev);CHKERRQ(ierr);
    //ierr = VecDestroy(&rhsVec);CHKERRQ(ierr);
    //ierr = MatDestroy(&mtx);CHKERRQ(ierr);
    //ierr = KSPDestroy(&ksp);CHKERRQ(ierr);

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

  ierr = MatAssemblyBegin(mtx,MAT_FINAL_ASSEMBLY);//CHKERRQ(ierr);
  ierr = MatAssemblyEnd(mtx,MAT_FINAL_ASSEMBLY);//CHKERRQ(ierr);
  
  //PetscViewerCreate(PETSC_COMM_WORLD, &viewer_matx);
  //PetscViewerDrawOpen();
  //PetscViewerSetFormat(viewer_matx, PETSC_VIEWER_ASCII_MATLAB);

  //MatView(mtx,PETSC_VIEWER_STDOUT_WORLD);

  VecZeroEntries(soln);

  VecAssemblyBegin(rhsVec);
  VecAssemblyEnd(rhsVec);

  VecAssemblyBegin(soln);
  VecAssemblyEnd(soln);

  computerTime.go(fct);

  //VecView(rhsVec, PETSC_VIEWER_STDOUT_WORLD);

  tstart = time(0);

  ierr = KSPSolve(ksp,rhsVec,soln);//CHKERRQ(ierr);

  ierr = KSPGetIterationNumber(ksp,&its);//CHKERRQ(ierr);

  KSPGetConvergedReason(ksp,&reason);

  if(reason<0)
  {
    printf("Divergence.\n");
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




int SolverPetsc::AssembleMatrixAndVector(vector<int>& row, vector<int>& col, MatrixXd& Klocal, VectorXd& Flocal)
{
  int ii, jj;
  for(ii=0;ii<row.size();ii++)
  {
    VecSetValues(rhsVec, 1, &row[ii], &Flocal(ii), ADD_VALUES);
    for(jj=0;jj<col.size();jj++)
    {
      //cout << ii << '\t' << jj << endl;
      MatSetValues(mtx, 1, &row[ii], 1, &col[jj], &(Klocal(ii,jj)), ADD_VALUES);
    }
  }

  return 1;
}



int SolverPetsc::AssembleMatrixAndVector(int start, int c1, vector<int>& vec1, vector<int>& vec2, MatrixXd& Klocal, VectorXd& Flocal)
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

/*
  if(STABILISED)
  {
  for(ii=0;ii<size2;ii++)
  {
    aa = start + vec2[ii];
    bb = size1 + ii;
    for(jj=0;jj<size2;jj++)
    {
      c = start+vec2[jj];
      MatSetValues(mtx, 1, &aa, 1, &c, &(Klocal(bb,size1+jj)), ADD_VALUES);
    }
  }
  }
*/
  return 0;
}




int SolverPetsc::AssembleMatrixAndVector(int start, int c1, vector<int>& vec1, MatrixXd& Klocal, VectorXd& Flocal)
{
  int ii, jj, size1;

  //printVector(vec1);
  //printVector(vec2);
  
  size1 = vec1.size();

  for(ii=0;ii<size1;ii++)
  {
    VecSetValue(rhsVec, vec1[ii], Flocal(ii), ADD_VALUES);
    //VecSetValues(rhsVec, 1, &vec1[ii], &Flocal(ii), ADD_VALUES);

    for(jj=0;jj<size1;jj++)
      MatSetValue(mtx, vec1[ii], vec1[jj], Klocal(ii,jj), ADD_VALUES);
      
      //MatSetValues(mtx, 1, &vec1[ii], 1, &vec1[jj], &(Klocal(ii,jj)), ADD_VALUES);
      //ierr = MatSetValue(mtx, aa, bb, Klocal(ii, jj), ADD_VALUES);
  }

  return 0;
}



int SolverPetsc::AssembleVector(int start, int c1, vector<int>& vec1, VectorXd& Flocal)
{
  int ii, jj;

  for(ii=0;ii<vec1.size();ii++)
  {
    VecSetValues(rhsVec, 1, &vec1[ii], &Flocal(ii), ADD_VALUES);
  }

  return 0;
}



int SolverPetsc::AssembleMatrixAndVectorCutFEM(int start, int c1, vector<int>& tempVec, vector<int>& forAssy, MatrixXd& Klocal, VectorXd& Flocal)
{
  int ii, jj, size1, r;

  //printVector(vec1);
  //printVector(vec2);

/*
  size1 = tempVec.size();

  for(ii=0;ii<size1;ii++)
  {
    r = start+forAssy[tempVec[ii]];

    rhsVec[r] += Flocal(ii);

    for(jj=0;jj<size1;jj++)
      mtx.coeffRef(r, start+forAssy[tempVec[jj]]) += Klocal(ii, jj);
  }
*/
//
  size1 = tempVec.size();

  for(ii=0;ii<size1;ii++)
  {
    r = forAssy[tempVec[ii]];

    VecSetValues(rhsVec, 1, &r, &Flocal(ii), ADD_VALUES);

    for(jj=0;jj<size1;jj++)
      MatSetValues(mtx, 1, &r, 1, &forAssy[tempVec[jj]], &(Klocal(ii,jj)), ADD_VALUES);
  }
//

  return 0;
}



int SolverPetsc::AssembleMatrixAndVectorCutFEM2(int start1, int start2, vector<int>& tempVec, 
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



int SolverPetsc::AssembleMatrixAndVectorCutFEM3(int start1, int start2, vector<int>& row, vector<int>& col, 
     vector<int>& forAssy1, vector<int>& forAssy2, MatrixXd& Klocal, VectorXd& Flocal1)
{
 
  return 0;
}




