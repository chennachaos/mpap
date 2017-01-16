
#ifndef  MY_POLYPRECOND_H
#define  MY_POLYPRECOND_H


#include <vector>
#include "../lib/DomainTree.h"

#include "headersEigen.h"

using std::vector;

using namespace std;
using namespace Eigen;
using namespace Eigen::internal;


namespace myIterSolvers {


typedef  int Index;
  
/*************************************************************

Polynominal preconditioner for solving the linear system
Ax = b

Input  --- Matrix A, rhs vector b
Output --- vector x

*************************************************************/


/** \internal Successive Overrelaxation (SOR) algorithm
  * \param mat The matrix A
  * \param rhs The right hand side vector b
  * \param x On input and initial solution, on output the computed solution.
  * \param iters On input the max number of iteration, on output the number of performed iterations.
  * \param tol_error On input the tolerance error, on output an estimation of the relative error.
  */



bool myPolyprecond_Dense(const MatrixXd& mat, const VectorXd& rhs, VectorXd& x, Index& Iters, 
           double& tol_error, int nTerms)
{
  using std::sqrt;
  using std::abs;

  typedef typename VectorXd::RealScalar RealScalar;
  typedef typename VectorXd::Scalar Scalar;

  RealScalar tol = tol_error, error_norm;

  Index maxIters = Iters;

  Index nRow = mat.cols();
  
  VectorXd y   = VectorXd::Zero(nRow);
  VectorXd z   = VectorXd::Zero(nRow);

  VectorXd  diagVec = mat.diagonal();
  VectorXd  diagInv = diagVec.cwiseInverse();
  
  //printVector(diagVec);
  //printVector(diagInv);

  VectorXd  resi   = rhs - mat*x;

  RealScalar rhs_norm = rhs.norm();

  if(rhs_norm == 0.0)
  {
    x.setZero();
    return true;
  }

  Index iter = 0, i, j;
  bool flag=false;
  //cout << " AAAAAAAAAAAAAAA " << endl;

/*
resi=b;
x=zeros(4,1);

for iter=1:10
    z = resi;
    for i=1:10
        y = resi;
        for j=1:i
            y = C*(invD*y) ;
        end
        z = z + y;
    end
    x = x + invD*z
    resi = b - A*x;
end
*/
  
  //resi = rhs;
  //x.setZero();

  //printVector(resi);

  while( iter < maxIters )
  {
    z = resi;
    for(i=0;i<nTerms;i++)
    {
      y = resi;
      for(j=0; j<=i; j++)
      {
        y = diagInv.cwiseProduct(y);
        y = diagVec.cwiseProduct(y) - mat*y;
        
        //printf(" y\n");        printVector(y);
      }
      //printf(" y\n");      printVector(y);

      z += y;
      //printf(" z\n");      printVector(z);
    }

    for( i=0; i<nRow; i++ )
      x[i] += diagInv[i]*z[i];
    
    //printVector(x);
    
    resi = rhs-mat*x;

    error_norm = resi.norm() / rhs_norm;
    cout << iter << '\t' << error_norm << endl;

    // Check for convergence.
    if( error_norm <= tol )
    {
      flag = 1;
      break;
    }
    ++iter;
  }

  Iters = iter;
  tol_error = error_norm;

  return flag;
}




bool myPolyprecond_Sparse(const SparseMatrixXd& mat, const VectorXd& rhs, VectorXd& x, Index& Iters, 
           double& tol_error, int&  nTerms)
{
  using std::sqrt;
  using std::abs;

  typedef typename VectorXd::RealScalar RealScalar;
  typedef typename VectorXd::Scalar Scalar;

  RealScalar tol = tol_error, error_norm;

  Index maxIters = Iters;

  Index nRow = mat.cols();
  
  VectorXd y   = VectorXd::Zero(nRow);
  VectorXd z   = VectorXd::Zero(nRow);

  VectorXd  diagVec = mat.diagonal();
  VectorXd  diagInv = diagVec.cwiseInverse();
  
  //printVector(diagVec);
  //printVector(diagInv);

  VectorXd  resi   = rhs - mat*x;

  RealScalar rhs_norm = rhs.norm();

  if(rhs_norm == 0.0)
  {
    x.setZero();
    return true;
  }

  Index iter = 0, i, j;
  bool flag=false;
  //cout << " AAAAAAAAAAAAAAA " << endl;


  //resi = rhs;
  //x.setZero();

  //printVector(resi);

  while( iter < maxIters )
  {
    z = resi;
    for(i=0;i<nTerms;i++)
    {
      y = resi;
      for(j=0; j<=i; j++)
      {
        y = diagInv.cwiseProduct(y);
        y = diagVec.cwiseProduct(y) - mat*y;
        
        //printf(" y\n");        printVector(y);
      }
      //printf(" y\n");      printVector(y);

      z += y;
      //printf(" z\n");      printVector(z);
    }

    for( i=0; i<nRow; i++ )
      x[i] += diagInv[i]*z[i];
    
    //printVector(x);
    
    resi = rhs-mat*x;

    error_norm = resi.norm() / rhs_norm;
    cout << iter << '\t' << error_norm << endl;

    // Check for convergence.
    if( error_norm <= tol )
    {
      flag = 1;
      break;
    }
    ++iter;
  }

  Iters = iter;
  tol_error = error_norm;

  return flag;
}






}

#endif


