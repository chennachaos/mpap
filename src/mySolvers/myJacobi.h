
#ifndef  MY_Jacobi_H
#define  MY_Jacobi_H


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

Jacobi iteration Method for solving the linear system

Ax = b

Input  --- Matrix A, rhs vector b
Output --- vector x

*************************************************************/


/** \internal Jacobi iteration algorithm
  * \param mat The matrix A
  * \param rhs The right hand side vector b
  * \param x On input and initial solution, on output the computed solution.
  * \param iters On input the max number of iteration, on output the number of performed iterations.
  * \param tol_error On input the tolerance error, on output an estimation of the relative error.
  */


template<typename MatrixType, typename VectorType>
bool myJacobi_Dense(const MatrixType& mat, const VectorType& rhs, VectorType& x, Index& Iters, 
           typename VectorType::RealScalar& tol_error, typename VectorType::RealScalar&  omega)
{
  using std::sqrt;
  using std::abs;

  typedef typename VectorType::RealScalar RealScalar;
  typedef typename VectorType::Scalar Scalar;

  RealScalar tol = tol_error, error_norm;

  Index maxIters = Iters;

  Index nRow = mat.cols();

  VectorType  xOld   = VectorType::Zero(nRow);

  RealScalar rhs_norm = rhs.norm();

  if(rhs_norm == 0.0)
  {
    x.setZero();
    return true;
  }

  Index iter = 0, i, j;
  bool flag=false;

  while( iter < maxIters )
  {
    xOld = x;

    for(i=0; i<nRow; i++ )
    {
      x[i] = rhs[i];

      for( j=0; j<i; j++ )
        x[i] -= mat(i,j) * xOld[j];

      for( j=i+1; j<nRow; j++ )
        x[i] -= mat(i,j) * xOld[j];

      x[i] /= mat(i,i);
    }

    error_norm = (rhs-mat*x).norm() / rhs_norm;
    //cout << iter << '\t' << error_norm << endl;

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



bool myJacobi_Sparse(const SparseMatrixXd& mat, const VectorXd& rhs, VectorXd& x, Index& Iters, 
          double& tol_error, double&  omega)
{
  using std::sqrt;
  using std::abs;

  //typedef typename VectorXd::RealScalar RealScalar;
  //typedef typename VectorXd::Scalar Scalar;

  typedef  double  RealScalar;
  typedef  double  Scalar;

  RealScalar tol = tol_error, error_norm;

  Index maxIters = Iters;

  Index nRow = mat.cols();

  VectorXd xOld   = VectorXd::Zero(nRow);

  VectorXd  diagVec = mat.diagonal();
  VectorXd  diagInv = diagVec.cwiseInverse();

  RealScalar rhs_norm = rhs.norm();

  if(rhs_norm == 0.0)
  {
    x.setZero();
    return true;
  }

  Index iter = 0, i, j;
  bool flag=false;
  //cout << " AAAAAAAAAAAAAAA " << endl;

  const Index* rowPtr     = mat.outerIndexPtr();
  const Index* colIndx    = mat.innerIndexPtr();
  const RealScalar* array = mat.valuePtr();
  
  //for(i=0;i<nRow;i++)
    //cout << i << '\t' << rhs[i] << endl;
/*
  while( iter < maxIters )
  {
    xOld = x;

    for(i=0; i<nRow; i++ )
    {
      x[i] = rhs[i];

      Index  k1  = rowPtr[i];
      Index  ikk = rowPtr[i+1]-rowPtr[i];
      Index  ik;
    
      //cout << k1 << '\t' << ikk << endl;
      for(ik=0; ik<ikk; ik++)
      {
        Index k2 = k1+ik;
        Index k = colIndx[k2];
        RealScalar  temp = array[k2];

        if( !(k == i) )
          x[i] -= temp * xOld[k];
      }

      x[i] *= diagInv[i];
    }

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
*/

  x = rhs;
  //for(i=0;i<2;i++)
    //x = diagInv.cwiseProduct(rhs);
    //x = diagInv.cwiseProduct(x);

  while( iter < maxIters )
  {
    //xOld = x;

    x += (omega*diagInv.cwiseProduct(x));

    error_norm = (rhs-mat*x).norm() / rhs_norm;
    //cout << iter << '\t' << error_norm << endl;

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


