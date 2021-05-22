#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#define EIGEN_SUPERLU_SUPPORT


#include "Debug.h"
#include "SolverEigen.h"
#include "SolverPetsc.h"
#include "SolverTime.h"
#include "ComputerTime.h"
#include "util.h"

//#include <Eigen/SuperLUSupport>
#include <Eigen/SparseExtra>
#include <Eigen/IterativeSolvers>



extern SolverTime      solverTime;
extern ComputerTime    computerTime;


using namespace std;
using namespace Eigen;



SolverEigen::SolverEigen()
{
  STABILISED = false;
  
  update_precond = 1;

  if (debug) cout << " SolverEigen constructor\n\n";
}


SolverEigen::~SolverEigen()
{
  if (debug)  cout << " SolverEigen destructor\n\n";

  free();
}


int SolverEigen::initialise(int p1, int p2, int p3)
{
   nRow = nCol = p3;

   soln.resize(nRow);
   soln.setZero();

   rhsVec   = soln;
   solnPrev = soln;

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

    SimplicialLDLT<SparseMatrix<double> > solver;

    //SuperLU<SparseMatrixXd > solver;
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

    //BiCGSTAB<SparseMatrixXd> solver;

    //GMRES<SparseMatrixXd, IncompleteLUT<double> > solver;

    //GMRES<SparseMatrixXd> solver;

    //ConjugateGradient<SparseMatrixXd, Lower, IncompleteLUT<double> >   solver;

    //solver.set_restart(100);
    //solver.setEigenv(10);

    solver.preconditioner().setDroptol(1.0e-3);
    solver.preconditioner().setFillfactor(1);

    solver.setMaxIterations(2000);
    solver.setTolerance(1.0e-8);

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

  solnPrev = soln;

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


