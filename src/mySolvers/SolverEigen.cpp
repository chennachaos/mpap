#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#define EIGEN_SUPERLU_SUPPORT


#include "Debug.h"
#include "SolverEigen.h"
#include "SolverPetsc.h"
#include "SolverTime.h"
#include "ComputerTime.h"
#include "util.h"

#include <Eigen/SuperLUSupport>
#include <Eigen/SparseExtra>
#include <Eigen/IterativeSolvers>



extern SolverTime      solverTime;
extern ComputerTime    computerTime;

/*
#define VIENNACL_WITH_EIGEN 1
#define VIENNACL_WITH_OPENMP 1

// ViennaCL headers
#include "viennacl/linalg/ilu.hpp"
#include "viennacl/linalg/cg.hpp"
#include "viennacl/linalg/bicgstab.hpp"
#include "viennacl/linalg/gmres.hpp"
#include "viennacl/linalg/jacobi_precond.hpp"
//#include "viennacl/linalg/spai.hpp"

#include "viennacl/io/matrix_market.hpp"
#include "viennacl/compressed_matrix.hpp"
//#include "Random.hpp"
//#include "vector-io.hpp"
typedef viennacl::compressed_matrix<double> SparseMatrixVCL;

//using namespace viennacl::linalg;
*/


using namespace std;
using namespace Eigen;



SolverEigen::SolverEigen()
{
  STABILISED = false;
  
  update_precond = 1;

  if (debug) cout << " SolverEigen constructor\n\n";

  //cout << " AAAAAAAAAA " << endl;
}


SolverEigen::~SolverEigen()
{
  if (debug)  cout << " SolverEigen destructor\n\n";

  free();
}


int SolverEigen::initialise(int p1, int p2, int p3)
{
   nRow = nCol = p3;

   //cout << " nRow = " << nRow << endl;

   soln.resize(nRow);
   soln.setZero();
   
   rhsVec   = soln;
   solnPrev = soln;

   //mtx.resize(nRow, nCol);

   //cout << " AAAAAAAAAA " << endl;
   //setSolverAndParameters();
   //cout << " AAAAAAAAAA " << endl;

  return 0;
}


int SolverEigen::setSolverAndParameters()
{
    ///////////////////////
    // Create the linear solver and set various options
    ///////////////////////


    ///////////////////////
    //Set operators. Here the matrix that defines the linear system
    //also serves as the preconditioning matrix.
    ///////////////////////

    return 0;
}


void SolverEigen::zeroMtx()
{
  //cout << " nRow = " << nRow << endl;
  //printVector(rhsVec);

  //cout << " nRow = " << nRow << endl;
  //cout << mtx << endl;
  //mtx.setZero();
  mtx *= 0.0;
  rhsVec.setZero();
  //cout << " nRow = " << nRow << endl;

  return;
}



int SolverEigen::free()
{
  return 0;
}


void SolverEigen::printInfo()
{
  //cout << "Eigen solver:  nRow = " << nRow << "\n";
  //cout << "               nnz  = " <<  << "\n\n"; 
  //printVector(rhsVec);
  //rhsVec.setZero();
  //cout << " nRow = " << nRow << endl;
  //cout << mtx << endl;

  return;
}


void SolverEigen::printMatrixPatternToFile()
{
    /*
    ofstream fout("matrix-pattern2.dat");

    if(fout.fail())
    {
      cout << " Could not open the Output file" << endl;
      exit(1);
    }

    fout << totalDOF << setw(10) << totalDOF << endl;

    for(int k=0; k<mtx.outerSize(); ++k)
    for(SparseMatrixXd::InnerIterator it(mtx,k); it; ++it)
      printf("%9d \t %9d \n", it.row(), it.col());

    //for(int ii=0;ii<nRow;ii++)
      //for(int jj=0;jj<nRow;jj++)
        //printf("%5d \t %5d \t %12.6f \n",mtx.coeffRef(ii,jj));

    fout.close();
    */

   FILE * pFile;
   int n;
   char name [100];

   cout << mtx.nonZeros() << endl;
   //pFile = fopen ("Stokes.dat","w");
   pFile = fopen ("Poisson.dat","w");

    fprintf(pFile, "%9d \t %9d \t %9d \n", nRow, nCol, 0 );

    for(int k=0; k<mtx.outerSize(); ++k)
    for(SparseMatrixXd::InnerIterator it(mtx,k); it; ++it)
      fprintf(pFile, "%9d \t %9d \t %20.16f \n", it.row(), it.col(), it.value());

   fclose (pFile);

   pFile = fopen ("rhsVec.dat","w");

    for(int k=0; k<nRow; ++k)
      fprintf(pFile, "%9d \t %14.8f \n", k, rhsVec(k));

   fclose (pFile);

   pFile = fopen ("refsoln.dat","w");

    for(int k=0; k<nRow; ++k)
      fprintf(pFile, "%9d \t %14.8f \n", k, soln(k));

   fclose (pFile);
  return;
}

  

void SolverEigen::printMatrix(int dig, int dig2, bool gfrmt, int indent, bool interactive)
{
  printInfo();
  
  cout << mtx << endl;
  printf("\n\n");

  return;
}

double SolverEigen::giveMatrixCoefficient(int row, int col)
{ 
  return  mtx.coeff(row,col);
}



int SolverEigen::factorise()
{
  char fct[] = "SolverEigen::factorise";

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


int  SolverEigen::solve()
{
  char fct[] = "SolverEigen::solve";

  time_t tstart, tend;

  if (currentStatus != FACTORISE_OK) { prgWarning(1,fct,"factorise matrix first!"); return 1; }
  

  //algoType = 1;
  
  //cout << mtx << endl;

//#ifdef EIGEN_PARALLELIZE
//#ifdef EIGEN_HAS_OPENMP
//printf("Eigen parallellize is on \n");
//#else
//printf("Eigen parallellize is off \n");
//#endif
  
  if(algoType == 1)
  {
    //cout << " Solving with Eigen::SimplicialLDLT " << endl;

    //SimplicialLDLT<SparseMatrix<double> > solver;

    SuperLU<SparseMatrixXd > solver;
    tstart = time(0);
  
    solver.compute(mtx);

    soln = solver.solve(rhsVec);

    //printVector(soln);
    //
    //computeConditionNumber();

    //VectorXd x0 = VectorXd::LinSpaced(nRow, 0.0, 1.0);

    //double  cond = myCondNum(mtx, x0, 50, solver);

    //printf("\n Matrix condition number = %12.6E \n\n\n", cond);

    //soln = solver.solveWithGuess(rhsVec, soln);
    //cout << solver.info() << '\t' << solver.error() << '\t' << solver.iterations() << endl;
    tend = time(0);
    //printf("It took %8.4f second(s) \n ", difftime(tend, tstart) );
  }
 
  if( algoType == 2 )
  {
    cout << " Solving with Eigen::BiCGSTAB " << endl;

    BiCGSTAB<SparseMatrixXd, IncompleteLUT<double> > solver;

    solver.preconditioner().setDroptol(1.0e-3);
    solver.preconditioner().setFillfactor(3);
    
    //BiCGSTAB<SparseMatrixXd> solver;

    //GMRES<SparseMatrixXd, IncompleteLUT<double> > solver;

    //GMRES<SparseMatrixXd> solver;

    //ConjugateGradient<SparseMatrixXd, Lower, IncompleteLUT<double> >   solver;

    
    //solver.set_restart(100);
    //solver.setEigenv(10);

    solver.setMaxIterations(500);
    solver.setTolerance(1.0e-10);

    tstart = time(0);
  
    //cout << " iiiiiiiiiiiiiii " << endl;
    solver.compute(mtx);
    //cout << " iiiiiiiiiiiiiii " << endl;
 
    soln = solver.solve(rhsVec);
  
    //printf("\n\n");
    //printVector(soln);

    //soln = solver.solveWithGuess(rhsVec, soln);

    tend = time(0); 
    printf("It took %8.4f second(s) \n ", difftime(tend, tstart) );

    cout << solver.info() << '\t' << solver.error() << '\t' << solver.iterations() << endl;
  }

  if( algoType == 3 )
  {
    /*
    cout << " Solving with ViennaCL::BiCGSTAB " << endl;
    
    solnPrev = soln;
    
    unsigned int aaaa = 8;

    viennacl::compressed_matrix<double, 1> vcl_compressed_matrix_1;
    //
    //SparseMatrixVCL     vcl_sparsematrix(nRow, nCol);
    viennacl::compressed_matrix<double>   vcl_sparsematrix(nRow, nCol, mtx.nonZeros());

    viennacl::vector<double>   vcl_rhs(nRow), vcl_soln(nRow), vcl_solnPrev(nRow);

    cout << " Copying the matrix and vector " << endl;

    viennacl::copy(mtx, vcl_sparsematrix);
    viennacl::copy(rhsVec, vcl_rhs);
    viennacl::copy(solnPrev, vcl_solnPrev);

    //
    cout << " Copied the matrix and vector. Solving the matrix system. " << endl;
    
    int  fillfactor=10;
    int  fill_in = (mtx.nonZeros()*fillfactor)/nRow + 1;
    cout << nRow << '\t' << fill_in << endl;
    if(fill_in > nRow)
      fill_in = nRow;

    // number of largest nonzero elements to keep in the L and the U part of the current row:
    int nnzL = fill_in*3;

    
    //cout << mtx << endl;
    //cout << "\n\n\n\n\n\n\n\n\n\n\n\n" << endl;
    //cout << vcl_sparsematrix << endl;
    
    tstart = time(0);

    // configuration of preconditioner:
    viennacl::linalg::chow_patel_tag  chow_patel_ilu_config;
    chow_patel_ilu_config.sweeps(2);       // three nonlinear sweeps
    chow_patel_ilu_config.jacobi_iters(50); // two Jacobi iterations per triangular 'solve' Rx=r
    // create and compute preconditioner:
    viennacl::linalg::chow_patel_ilu_precond<viennacl::compressed_matrix<double> > chow_patel_ilu(vcl_sparsematrix, chow_patel_ilu_config);


    viennacl::linalg::ilut_tag  ilut_config(nnzL, 1.0e-3);
    viennacl::linalg::ilu0_tag  ilu0_config;
    viennacl::linalg::jacobi_tag  jacobi_config;
    //viennacl::linalg::spai_tag spai_config(1.0e-3, 3, 5.0e-2);

    viennacl::linalg::ilu0_precond<viennacl::compressed_matrix<double> >  vcl_ilu0(vcl_sparsematrix, ilu0_config);
    //viennacl::linalg::ilut_precond<viennacl::compressed_matrix<double> >  vcl_ilut(vcl_sparsematrix, ilut_config);
    //viennacl::linalg::jacobi_precond<SparseMatrixVCL >  vcl_jacobi(vcl_sparsematrix, jacobi_config);

    // compute block-ILU preconditioner using ILU0 for each block:
    viennacl::linalg::block_ilu_precond<viennacl::compressed_matrix<double>, viennacl::linalg::ilu0_tag> vcl_block_ilu0(vcl_sparsematrix, ilu0_config, 1);

    // compute block-ILU preconditioner using ILUT for each block:
    //viennacl::linalg::block_ilu_precond<viennacl::compressed_matrix<double>, viennacl::linalg::ilut_tag> vcl_block_ilut(vcl_sparsematrix, ilut_config);

    //viennacl::linalg::spai_precond<viennacl::compressed_matrix<double> > vcl_spai(vcl_sparsematrix, spai_config);

    //cout << " jjjjjjjjjjjjj " << endl;

    //viennacl::linalg::cg_tag  custom_tag(1.0e-12, 1000);
    viennacl::linalg::bicgstab_tag  custom_tag(1e-10, 500);
    //viennacl::linalg::gmres_tag  custom_tag(1.0e-10, 500, 50);

    //viennacl::linalg::bicgstab_solver<viennacl::vector<double> > my_bicgstab_solver(custom_tag);
    //my_bicgstab_solver.set_initial_guess(vcl_solnPrev);

    //tstart = time(0);

    //vcl_soln = viennacl::linalg::solve(vcl_sparsematrix, vcl_rhs, custom_tag, chow_patel_ilu); // preconditioner here

    //vcl_soln = viennacl::linalg::solve(vcl_sparsematrix, vcl_rhs, custom_tag, vcl_ilu0);
    //vcl_soln = viennacl::linalg::solve(vcl_sparsematrix, vcl_rhs, custom_tag, vcl_ilut);
    //vcl_soln = viennacl::linalg::solve(vcl_sparsematrix, vcl_rhs, custom_tag, vcl_jacobi);

    vcl_soln = viennacl::linalg::solve(vcl_sparsematrix, vcl_rhs, custom_tag, vcl_block_ilu0);
    //vcl_soln = viennacl::linalg::solve(vcl_sparsematrix, vcl_rhs, custom_tag, vcl_block_ilut);
    
    //vcl_soln = viennacl::linalg::solve(vcl_sparsematrix, vcl_rhs, custom_tag, vcl_spai);

    //vcl_soln = my_bicgstab_solver(vcl_sparsematrix, vcl_rhs, vcl_ilu0);

 
    tend = time(0); 
    printf("It took %8.4f second(s) \n ", difftime(tend, tstart) );

    cout << "No. of iters: " << custom_tag.iters() << endl;
    cout << "Est. error: " << custom_tag.error() << endl;

    //tstart = time(0);

    cout << " jjjjjjjjjjjjj " << endl;

    viennacl::copy(vcl_soln, soln);
    */
  }

  if( algoType == 4 )
  {
    //solverSchurCG();
    solverSchurBiCGSTAB();
    //solverSchurGMRES();
  }

  if( algoType == 5 )
  {
    solverUzawaType1();
    //solverSchurGMRES();
  }

  if( algoType == 6 )
  {
    SolverPetsc  solverPetsc;
    cout << " Petsc solver " << endl;

    solverPetsc.solveSerial(mtx, rhsVec, soln);
    //solverPetsc.solveParallel(mtx, rhsVec, soln);
    //printMatrixPatternToFile();
  }

  int niter=2000;
  double  tol=1.0e-10;
  bool info=false;


  solnPrev = soln;

  //printf("\n\n");
  //printVector(soln);

  //for(int ii=0;ii<nRow;ii++)
    //RHS[ii] = soln(ii);
  
  //cout << "SolverEigen::solve()  took "<< difftime(tend, tstart) <<" second(s)."<< endl;

  //solverTime.total -= solverTime.solve;
  //solverTime.solve += computerTime.stop(fct);
  //solverTime.total += solverTime.solve;

  //if(INFO[0] != 0) return NULL;
  
  if(checkIO)
  {
    // search for "nan" entries in solution vector

    //if (prgNAN(RHS,N)) prgError(1,fct,"nan entry in solution vector!");
  }

  return 0;
}






int SolverEigen::factoriseAndSolve()
{
  if(currentStatus != ASSEMBLY_OK)
  {
    cerr << " assemble matrix first! " << endl;
    return 1;
  }
  
  factorise();

  return solve();
}



int SolverEigen::assembleMatrixAndVector(vector<int>& row, vector<int>& col, MatrixXd& Klocal, VectorXd& Flocal)
{
  int ii, jj;
  for(ii=0;ii<row.size();ii++)
  {
    rhsVec[row[ii]] += Flocal(ii);
    for(jj=0;jj<col.size();jj++)
    {
      mtx.coeffRef(row[ii], col[jj]) += Klocal(ii,jj);
    }
  }

  return 0;
}




int SolverEigen::assembleMatrixAndVector(int start, int c1, vector<int>& vec1, vector<int>& vec2, MatrixXd& Klocal, VectorXd& Flocal)
{
  // subroutine for mixed formulation
  // vec1 - for variable 'u'
  // vec2 - for variable 'p'

  int ii, jj, aa, bb, size1, size2;

  //printVector(vec1);
  //printVector(vec2);
  
  size1 = vec1.size();
  size2 = vec2.size();

  for(ii=0;ii<size1;ii++)
  {
    rhsVec[vec1[ii]] += Flocal(ii);
    for(jj=0;jj<size1;jj++)
       mtx.coeffRef(vec1[ii], vec1[jj]) += Klocal(ii, jj);

    for(jj=0;jj<size2;jj++)
    {
       aa = start + vec2[jj];
       bb = size1 + jj;
       mtx.coeffRef(vec1[ii], aa) += Klocal(ii, bb);
       mtx.coeffRef(aa, vec1[ii]) += Klocal(bb, ii);
    }
  }

  for(ii=0;ii<size2;ii++)
  {
    aa = start + vec2[ii];
    rhsVec[aa] += Flocal(size1+ii);
  }

  if(STABILISED)
  {
    for(ii=0;ii<size2;ii++)
    {
      aa = start + vec2[ii];
      bb = size1 + ii;
      for(jj=0;jj<size2;jj++)
      {
        mtx.coeffRef(aa, start+vec2[jj]) += Klocal(bb, size1+jj);
      }
    }
  }

  return 0;
}




int SolverEigen::assembleMatrixAndVector(int start, int c1, vector<int>& forAssy, MatrixXd& Klocal, VectorXd& Flocal)
{
  int ii, jj, aa, bb, size1, r, c;

  size1 = forAssy.size();

  /*
  for(ii=0;ii<size1;ii++)
  {
    rhsVec[forAssy[ii]] += Flocal(ii);

    for(jj=0;jj<size1;jj++)
      mtx.coeffRef(forAssy[ii], forAssy[jj]) += Klocal(ii, jj);
  }
  */

  for(ii=0; ii<size1; ii++)
  {
    aa = forAssy[ii];
    if( aa != -1 )
    {
      r = start + aa;
      rhsVec[r] += Flocal(ii);

      for(jj=0; jj<size1; jj++)
      {
        bb = forAssy[jj];
        if( bb != -1 )
          mtx.coeffRef(r, start+bb) += Klocal(ii,jj);
      }
    }
  }

  return 0;
}



int SolverEigen::assembleVector(int start, int c1, vector<int>& vec1, VectorXd& Flocal)
{
  int ii, jj;

  for(ii=0;ii<vec1.size();ii++)
  {
    rhsVec[vec1[ii]] += Flocal(ii);
  }

  return 0;
}



int SolverEigen::assembleMatrixAndVectorCutFEM(int start, int c1, vector<int>& grid2cutfem_DOF, vector<int>& forAssy, MatrixXd& Klocal, VectorXd& Flocal)
{
  int ii, jj, size1, r;

  //printVector(vec1);
  //printVector(vec2);

/*
  size1 = grid2cutfem_DOF.size();

  for(ii=0;ii<size1;ii++)
  {
    r = start+forAssy[grid2cutfem_DOF[ii]];

    rhsVec[r] += Flocal(ii);

    for(jj=0;jj<size1;jj++)
      mtx.coeffRef(r, start+forAssy[grid2cutfem_DOF[jj]]) += Klocal(ii, jj);
  }
*/
//
  size1 = grid2cutfem_DOF.size();

  for(ii=0;ii<size1;ii++)
  {
    r = forAssy[grid2cutfem_DOF[ii]];

    rhsVec[r] += Flocal(ii);

    for(jj=0;jj<size1;jj++)
      mtx.coeffRef(r, forAssy[grid2cutfem_DOF[jj]]) += Klocal(ii, jj);
  }
//

/*
  size1 = grid2cutfem_DOF.size();

  for(ii=0;ii<size1;ii++)
  {
    r = forAssy[ii];

    rhsVec[r] += Flocal(ii);

    for(jj=0;jj<size1;jj++)
      mtx.coeffRef(r, forAssy[jj]) += Klocal(ii, jj);
  }
*/

  return 0;
}



int SolverEigen::assembleMatrixAndVectorCutFEM2(int start1, int start2, vector<int>& grid2cutfem_DOF, 
     vector<int>& forAssy1, vector<int>& forAssy2, MatrixXd& Klocal, VectorXd& Flocal1, VectorXd& Flocal2)
{
  // to assemble coupling matrices arriving from 
  // jump conditions in the gradients across an interface

  int ii, jj, size1, kk, r, c;

  //printVector(vec1);
  //printVector(vec2);
  
  size1 = grid2cutfem_DOF.size();

  for(ii=0;ii<size1;ii++)
  {
    // force vector for domain #2
    r = start2 + forAssy2[grid2cutfem_DOF[ii]];

    rhsVec[r] += Flocal2(ii);

    // force vector for domain #1
    r = start1 + forAssy1[grid2cutfem_DOF[ii]];

    rhsVec[r] += Flocal1(ii);

    // coupling matrix
    for(jj=0;jj<size1;jj++)
    {
      c = start2 + forAssy2[grid2cutfem_DOF[jj]];
      mtx.coeffRef(r, c) += Klocal(ii, jj);
      mtx.coeffRef(c, r) += Klocal(ii, jj);
    }
  }

  return 0;
}



int SolverEigen::assembleMatrixAndVectorCutFEM3(int start1, int start2, vector<int>& row, vector<int>& col, 
     vector<int>& forAssy1, vector<int>& forAssy2, MatrixXd& Klocal, VectorXd& Flocal1)
{
 
  return 0;
}


