
#ifndef  MY_SOR_H
#define  MY_SOR_H


#include <vector>

#include "headersEigen.h"

using std::vector;

using namespace std;
using namespace Eigen;
using namespace Eigen::internal;


namespace myIterSolvers {


typedef  int Index;
  
/*************************************************************

Successive Overrelaxation (SOR) Method for solving the linear system

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


template<typename MatrixType, typename VectorType>
bool mySOR_Dense(const MatrixType& mat, const VectorType& rhs, VectorType& x, Index& Iters, 
           typename VectorType::RealScalar& tol_error, typename VectorType::RealScalar&  omega)
{
  using std::sqrt;
  using std::abs;

  typedef typename VectorType::RealScalar RealScalar;
  typedef typename VectorType::Scalar Scalar;

  RealScalar tol = tol_error, error_norm;

  Index maxIters = Iters;

  Index nRow = mat.cols();

  VectorType xOld   = VectorType::Zero(nRow);
  
  RealScalar rhs_norm = rhs.norm();

  if(rhs_norm == 0.0)
  {
    x.setZero();
    return true;
  }

  Index iter = 0, i, j;
  bool flag=false;
  //cout << " AAAAAAAAAAAAAAA " << endl;

  while( iter < maxIters )
  {
    xOld = x;

    //  Do the Gauss-Seidel computation.
    //
    for(i=0; i<nRow; i++ )
    {
      x[i] = rhs[i];

      for( j=0; j<i; j++ )
      {
        x[i] -= mat(i,j) * x[j];
      }

      for( j=i+1; j<nRow; j++ )
      {
        x[i] -= mat(i,j) * xOld[j];
      }

      x[i] /= mat(i,i);
    }

    //  Use W to blend the Gauss-Seidel update with the old solution.
    //
    //for( i=0; i<nRow; i++ )
      //x[i] = (1.0 - omega) * xOld[i] + omega*x[i];

    error_norm = (rhs-mat*x).norm() / rhs_norm;
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

  return flag;
}



bool mySOR_Sparse(const SparseMatrixXd& mat, const VectorXd& rhs, VectorXd& x, Index& Iters, 
          double& tol_error, double&  omega)
{
  using std::sqrt;
  using std::abs;

  typedef typename VectorXd::RealScalar RealScalar;
  typedef typename VectorXd::Scalar Scalar;

  RealScalar tol = tol_error, error_norm;

  Index maxIters = Iters;

  Index nRow = mat.cols();

  VectorXd  xOld   = VectorXd::Zero(nRow);
  VectorXd  diagInv = mat.diagonal().cwiseInverse();

  VectorXd  resi   = rhs;

  RealScalar rhs_norm = rhs.norm();

  if(rhs_norm == 0.0)
  {
    x.setZero();
    return true;
  }

  Index iter = 0, i, j;
  bool flag=false;
  //cout << " AAAAAAAAAAAAAAA " << endl;
  //cout << rhs_norm << endl;

  const Index* rowPtr     = mat.outerIndexPtr();
  const Index* colIndx    = mat.innerIndexPtr();
  const RealScalar* array = mat.valuePtr();
  
  Index k1, ik, ikk, k2, k;
  RealScalar  temp;

  while( iter < maxIters )
  {
    xOld = x;

    //  Do the Gauss-Seidel computation.
    //
    for(i=0; i<nRow; i++ )
    {
      x[i] = resi[i];

      k1 = rowPtr[i];
      ikk = (rowPtr[i+1]-rowPtr[i]);

      //cout << k1 << '\t' << ikk << endl;
      for(ik=0; ik<ikk; ik++)
      {
        k2 = k1+ik;
        k = colIndx[k2];
        temp = array[k2];

        if( k<i )
          x[i] -= temp * x[k];
        else if( k>i )
          x[i] -= temp * xOld[k];
      }

      x[i] *= diagInv[i];
      //x[i] /= mat.coeffRef(i,i);
      //cout << i << '\t' << mat.coeffRef(i,i) << endl;
    }

    //  Use W to blend the Gauss-Seidel update with the old solution.
    //
    // here Gauss-Siedel method is used
    for( i=0; i<nRow; i++ )
      x[i] = (1.0 - omega) * xOld[i] + omega*x[i];

    resi = rhs - mat*x;
    
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


