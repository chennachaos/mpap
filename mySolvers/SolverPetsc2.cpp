
#include "Debug.h"
#include "SolverPetsc.h"
#include "petscmat.h"
#include "petscksp.h"
#include "SolverTime.h"
#include "ComputerTime.h"


extern SolverTime      solverTime;
extern ComputerTime    computerTime;



int SolverPetsc::solveSerial(SparseMatrixXd& matEigen, VectorXd& rhsEigen, VectorXd& solnEigen)
{
    nRow = matEigen.rows();
    nCol = nRow;
   
    cout << " nRow = " << nRow << endl;

    time_t tstart, tend; 

    PetscInt its;

    vector<PetscInt> ix(nRow), nnzVec(nRow);

    KSPConvergedReason reason;
    PetscInt ii,  jj;
    PetscScalar  data;
    PetscScalar   zero = 0.0;

    //////////////////////////////////////////////////////
    //
    // Initialize Petsc
    //
    //////////////////////////////////////////////////////

    cout << " AAAAAAAAAA " << endl;

    PetscInitialize(NULL,NULL,(char *)0,NULL);
    //PetscInitialize(NULL, NULL, "petsc_options.dat", NULL);

    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    //if(size != 1)
      //SETERRQ(1,"This is a uniprocessor example only!");

    //ierr = PetscOptionsGetInt(PETSC_NULL, "-n", &nRow, PETSC_NULL);CHKERRQ(ierr);

    cout << " AAAAAAAAAA " << endl;

    VecCreate(PETSC_COMM_WORLD, &soln);
    //VecCreate(PETSC_COMM_WORLD, &solnPrev);
    //VecCreate(PETSC_COMM_WORLD, &rhsVec);
    //ierr = PetscObjectSetName((PetscObject) x, "Solution");CHKERRQ(ierr);

    cout << " BBBBBBBBBB " << endl;

    ierr = VecSetSizes(soln, PETSC_DECIDE, nRow);CHKERRQ(ierr);
    //ierr = VecSetSizes(solnPrev, PETSC_DECIDE, nRow);CHKERRQ(ierr);
    //ierr = VecSetSizes(rhsVec, PETSC_DECIDE, nRow);CHKERRQ(ierr);
    cout << " AAAAAAAAAA " << endl;
    ierr = VecSetFromOptions(soln);CHKERRQ(ierr);
    ierr = VecDuplicate(soln, &rhsVec);CHKERRQ(ierr);
    ierr = VecDuplicate(soln, &solnPrev);CHKERRQ(ierr);
    
    //////////////////////////////////////////////////////
    //
    // create and setup the matrix
    //
    //////////////////////////////////////////////////////

    cout << " Creating and setting Petsc matrix " << endl;

    //MatCreate(PETSC_COMM_WORLD, &mtx);

    //ierr = MatSetSizes(mtx,PETSC_DECIDE,PETSC_DECIDE,nRow,nRow);CHKERRQ(ierr);

    for(ii=0;ii<nRow;ii++)
      nnzVec[ii] = nCol/20;

    //ii = matEigen.nonZeros();
    ii = 2;
    MatCreateSeqAIJ(PETSC_COMM_WORLD, nRow, nCol, ii, &nnzVec[0], &mtx);
    //MatCreateSeqAIJ(PETSC_COMM_WORLD, nRow, nRow, ii, NULL, &mtx);

    //ierr = MatSetFromOptions(mtx);CHKERRQ(ierr);
    
    cout << " AAAAAAAAAA " << endl;
    //cout << matEigen << endl;
    cout << " AAAAAAAAAA " << endl;

    MatZeroEntries(mtx);
    
    matEigen.makeCompressed();

    for(int k=0; k<matEigen.outerSize(); ++k)
    {
      for(SparseMatrixXd::InnerIterator it(matEigen,k); it; ++it)
      {
        ii = it.row();
        jj = it.col();
        data = it.value();
        
        //cout << it.row() << '\t' << it.col() << '\t' << it.value() << endl;
        //cout << ii << '\t' << jj << '\t' << data << endl;

        MatSetValue(mtx, ii, jj, data, ADD_VALUES);
      }
      //cout << endl;      cout << endl;
    }

    ierr = MatAssemblyBegin(mtx,MAT_FINAL_ASSEMBLY);//CHKERRQ(ierr);
    ierr = MatAssemblyEnd(mtx,MAT_FINAL_ASSEMBLY);//CHKERRQ(ierr);

    cout << " Creating and setting Petsc matrix DONE " << endl;
    
    //VecSetValues(rhsVec, nRow, rhsEigen(ii), INSERT_VALUES);

    //VecSet(rhsVec, zero);
    for(ii=0;ii<nRow;ii++)
    {
      //cout << ii << endl;
      VecSetValue(rhsVec, ii, rhsEigen(ii), INSERT_VALUES);
      VecSetValue(soln, ii, 0.0, INSERT_VALUES);
      ix[ii] = ii;
    }

    cout << " AAAAAAAAAA " << endl;

    VecAssemblyBegin(rhsVec);
    VecAssemblyEnd(rhsVec);

    VecAssemblyBegin(soln);
    VecAssemblyEnd(soln);

    cout << " Creating and setting KSP " << endl;

    //////////////////////////////////////////////////////
    //
    // Create the linear solver and set various options
    //
    //////////////////////////////////////////////////////

    ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);CHKERRQ(ierr);

    ///////////////////////
    //Set operators. Here the matrix that defines the linear system
    //also serves as the preconditioning matrix.
    ///////////////////////

    ierr = KSPSetOperators(ksp,mtx,mtx);CHKERRQ(ierr);

    //ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

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
    //KSPSetType(ksp, KSPCG);
    //KSPSetType(ksp, KSPGMRES);
    KSPSetType(ksp, KSPBCGS);
    //KSPSetType(ksp, KSPBICG);
    //KSPSetType(ksp, KSPLSQR);
    //KSPGMRESSetRestart(ksp, 100);

    ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);

    //ierr = PCSetType(pc,PCJACOBI);CHKERRQ(ierr);
    ierr = PCSetType(pc,PCILU);CHKERRQ(ierr);
    //ierr = PCSetType(pc,PCLU);CHKERRQ(ierr);
    //ierr = PCSetType(pc,PCCholesky);CHKERRQ(ierr);
    //ierr = PCSetType(pc,PCNONE);CHKERRQ(ierr);
    //
    //PCFactorSetFill(pc, 8);
    //PCFactorsetLevels(pc, 5);

    ierr = KSPSetTolerances(ksp,1.0e-16,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);

    /*
    Set runtime options, e.g.,
    -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
    These options will override those specified above as long as
    KSPSetFromOptions() is called _after_ any other customization
    routines.
    */

    //ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

    //KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
    //KSPSetInitialGuessNonzero(ksp, PETSC_FALSE);

    //PetscViewerCreate(PETSC_COMM_WORLD, &viewer_matx);
    //PetscViewerDrawOpen();
    //PetscViewerSetFormat(viewer_matx, PETSC_VIEWER_ASCII_MATLAB);

    //MatView(mtx,PETSC_VIEWER_STDOUT_WORLD);

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
      printf("Convergence in %d iterations.\n",(int)its);
    }

    ierr = PetscPrintf(PETSC_COMM_WORLD,"Iterations %5D\n\n\n", its);//CHKERRQ(ierr);

    //ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);//CHKERRQ(ierr);

    //VecView(soln, PETSC_VIEWER_STDOUT_WORLD);

    VecGetValues(soln, nRow, &ix[0], &solnEigen(0));

    //for(ii=0;ii<nRow;ii++)
      //solnEigen[ii] =  soln[ii] ;

    tend = time(0); 
  
    cout << "SolverPetsc::solve()  took "<< difftime(tend, tstart) <<" second(s)."<< endl;

    //solverTime.total -= solverTime.solve;
    //solverTime.solve += computerTime.stop(fct);
    //solverTime.total += solverTime.solve;

//  ierr = PetscFinalize();CHKERRQ(ierr);

  return 0;
}



int SolverPetsc::solveParallel(SparseMatrixXd& matEigen, VectorXd& rhsEigen, VectorXd& solnEigen)
{
    nRow = matEigen.rows();
    nCol = nRow;
   
    cout << " nRow = " << nRow << endl;

    time_t tstart, tend; 

    vector<PetscInt>  ix(nRow), nnzVec(nRow);
    PetscInt its;
    KSPConvergedReason reason;
    PetscInt ii,  jj;
    PetscScalar  data;
    PetscScalar   zero = 0.0;

    //////////////////////////////////////////////////////
    //
    // Initialize Petsc
    //
    //////////////////////////////////////////////////////

    cout << " AAAAAAAAAA " << endl;

    //PetscInitialize(NULL,NULL,(char *)0,NULL);
    //PetscInitialize(NULL, NULL, "petsc_options.dat", NULL);

    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    cout << " AAAAAAAAAA " << endl;

    VecCreate(PETSC_COMM_WORLD, &soln);
    //VecCreate(PETSC_COMM_WORLD, &solnPrev);
    //VecCreate(PETSC_COMM_WORLD, &rhsVec);
    //ierr = PetscObjectSetName((PetscObject) x, "Solution");CHKERRQ(ierr);

    cout << " BBBBBBBBBB " << endl;

    ierr = VecSetSizes(soln, PETSC_DECIDE, nRow);CHKERRQ(ierr);
    //ierr = VecSetSizes(solnPrev, PETSC_DECIDE, nRow);CHKERRQ(ierr);
    //ierr = VecSetSizes(rhsVec, PETSC_DECIDE, nRow);CHKERRQ(ierr);
    cout << " AAAAAAAAAA " << endl;
    ierr = VecSetFromOptions(soln);CHKERRQ(ierr);
    ierr = VecDuplicate(soln, &rhsVec);CHKERRQ(ierr);
    ierr = VecDuplicate(soln, &solnPrev);CHKERRQ(ierr);
    
    //////////////////////////////////////////////////////
    //
    // create and setup the matrix
    //
    //////////////////////////////////////////////////////

    cout << " Creating and setting Petsc matrix " << endl;

    MatCreate(PETSC_COMM_WORLD, &mtx);

    ierr = MatSetSizes(mtx,PETSC_DECIDE,PETSC_DECIDE,nRow,nRow);CHKERRQ(ierr);

    for(ii=0;ii<nRow;ii++)
      nnzVec[ii] = nCol-1;


    //ii = matEigen.nonZeros();
    ii = 2;
    //MatCreateSeqAIJ(PETSC_COMM_WORLD, nRow, nCol, ii, nnzVec, &mtx);
    //MatCreateSeqAIJ(PETSC_COMM_WORLD, nRow, nRow, ii, NULL, &mtx);

    //ierr = MatSetFromOptions(mtx);CHKERRQ(ierr);

    ierr = MatMPIAIJSetPreallocation(mtx,5,PETSC_NULL,5,PETSC_NULL);CHKERRQ(ierr);
    
    cout << " AAAAAAAAAA " << endl;
    //cout << matEigen << endl;
    cout << " AAAAAAAAAA " << endl;

    MatZeroEntries(mtx);
    
    matEigen.makeCompressed();

    for(int k=0; k<matEigen.outerSize(); ++k)
    {
      for(SparseMatrixXd::InnerIterator it(matEigen,k); it; ++it)
      {
        ii = it.row();
        jj = it.col();
        data = it.value();
        
        //cout << it.row() << '\t' << it.col() << '\t' << it.value() << endl;
        //cout << ii << '\t' << jj << '\t' << data << endl;

        MatSetValue(mtx, ii, jj, data, ADD_VALUES);
      }
      //cout << endl;      cout << endl;
    }

    ierr = MatAssemblyBegin(mtx,MAT_FINAL_ASSEMBLY);//CHKERRQ(ierr);
    ierr = MatAssemblyEnd(mtx,MAT_FINAL_ASSEMBLY);//CHKERRQ(ierr);

    cout << " Creating and setting Petsc matrix DONE " << endl;
    
    //VecSetValues(rhsVec, nRow, rhsEigen(ii), INSERT_VALUES);

    //VecSet(rhsVec, zero);
    for(ii=0;ii<nRow;ii++)
    {
      //cout << ii << endl;
      VecSetValue(rhsVec, ii, rhsEigen(ii), INSERT_VALUES);
      VecSetValue(soln, ii, 0.0, INSERT_VALUES);
      ix[ii] = ii;
    }

    cout << " AAAAAAAAAA " << endl;

    VecAssemblyBegin(rhsVec);
    VecAssemblyEnd(rhsVec);

    VecAssemblyBegin(soln);
    VecAssemblyEnd(soln);

    cout << " Creating and setting KSP " << endl;

    //////////////////////////////////////////////////////
    //
    // Create the linear solver and set various options
    //
    //////////////////////////////////////////////////////

    ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);CHKERRQ(ierr);

    ///////////////////////
    //Set operators. Here the matrix that defines the linear system
    //also serves as the preconditioning matrix.
    ///////////////////////

    ierr = KSPSetOperators(ksp, mtx, mtx);CHKERRQ(ierr);

    //ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

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
    //KSPSetType(ksp, KSPCG);
    //KSPSetType(ksp, KSPGMRES);
    KSPSetType(ksp, KSPBCGS);
    //KSPSetType(ksp, KSPBICG);
    //KSPSetType(ksp, KSPLSQR);
    //KSPGMRESSetRestart(ksp, 100);

    ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);

    //ierr = PCSetType(pc,PCJACOBI);CHKERRQ(ierr);
    ierr = PCSetType(pc,PCILU);CHKERRQ(ierr);
    //ierr = PCSetType(pc,PCLU);CHKERRQ(ierr);
    //ierr = PCSetType(pc,PCCholesky);CHKERRQ(ierr);
    //ierr = PCSetType(pc,PCNONE);CHKERRQ(ierr);
    //
    //PCFactorSetFill(pc, 8);
    //PCFactorsetLevels(pc, 5);

    ierr = KSPSetTolerances(ksp,1.0e-16,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);

    /*
    Set runtime options, e.g.,
    -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
    These options will override those specified above as long as
    KSPSetFromOptions() is called _after_ any other customization
    routines.
    */

    //ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

    //KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
    //KSPSetInitialGuessNonzero(ksp, PETSC_FALSE);

    //PetscViewerCreate(PETSC_COMM_WORLD, &viewer_matx);
    //PetscViewerDrawOpen();
    //PetscViewerSetFormat(viewer_matx, PETSC_VIEWER_ASCII_MATLAB);

    //MatView(mtx,PETSC_VIEWER_STDOUT_WORLD);

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
      printf("Convergence in %d iterations.\n",(int)its);
    }

    ierr = PetscPrintf(PETSC_COMM_WORLD,"Iterations %5D\n\n\n", its);//CHKERRQ(ierr);

    //ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);//CHKERRQ(ierr);

    //VecView(soln, PETSC_VIEWER_STDOUT_WORLD);

    VecGetValues(soln, nRow, &ix[0], &solnEigen(0));

    //for(ii=0;ii<nRow;ii++)
      //solnEigen[ii] =  soln[ii] ;

    tend = time(0); 
  
    cout << "SolverPetsc::solve()  took "<< difftime(tend, tstart) <<" second(s)."<< endl;

    //solverTime.total -= solverTime.solve;
    //solverTime.solve += computerTime.stop(fct);
    //solverTime.total += solverTime.solve;

//  ierr = PetscFinalize();CHKERRQ(ierr);

  return 0;
}



