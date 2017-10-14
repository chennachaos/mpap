
#ifndef MY_COND_NUMBER_H
#define MY_COND_NUMBER_H


#include "headersBasic.h"

#include <Eigen/SuperLUSupport>

typedef  SparseMatrix<double, RowMajor>  SparseMatrixXd;


template <typename VectorType>
double  infNormVector(const VectorType& x)
{
  return  x.cwiseAbs().maxCoeff();
}



//!  This function is used to find the condition number of a square matrix
/*
*
*/
//template <typename MatrixType>

double  myCondNum(const SparseMatrixXd& A, const VectorXd& x0, int nIter, SuperLU<SparseMatrixXd>& solver)
{
  int i,j,n = A.rows();

  VectorXd  x(n), xNew(n);

  double  lamOld=0.0, lamMax=0.0, lamMin=0.0, cond, tol=1.0e-4;


    // compute the largest eigenvalue of A

    x = x0;

    for(i=0; i<nIter; i++)
    {
      xNew = A*x;

      // using power iteration method
      //lam1 = xNew.cwiseAbs().maxCoeff()/x.cwiseAbs().maxCoeff();
      lamMax = infNormVector(xNew)/infNormVector(x);

      // using Rayleigh quotient method

      //lam = xNew.dot(xNew)/x.norm();

      printf(" \t %5d \t %12.8f \n ", i, lamMax);

      if(abs(lamMax-lamOld) < tol)
        break;

      lamOld = lamMax;

      x = xNew;
    }

    //printf(" Largest eigenvalue = %12.6f \n", lam1);

    // compute the smallest eigenvalue of A
    // which is nothing but the largest eigenvalue of inv(A)
    // using power iteration method

    //SuperLU<SparseMatrix<double> > solver;
    
    //SuperLU<SparseMatrixXd> solver;
    
    //SimplicialLDLT<SparseMatrix<double> > solver;
    
    //FullPivLU<MatrixXd> solver(A);
  
    //solver.compute(A);

    x = x0.cwiseInverse();
    x = x0;

    lamMin = 0.0;
    for(i=0; i<nIter; i++)
    {
      xNew = solver.solve(x);

      //lamMin = infNormVector(x1)/infNormVector(x);
      lamMin = xNew.cwiseAbs().maxCoeff()/x.cwiseAbs().maxCoeff();

      printf(" \t %5d \t %12.8f \n ", i, lamMin);

      //if( (abs(lam2-tmp) < 1.0e-3) || (abs(lam2) > 1.0e8) )
      if( (abs(lamOld-lamMin) < tol) )
        break;

      lamOld = lamMin;

      x = xNew;
    }
    lamMin = 1.0/lamMin;

    cond = abs(lamMax/lamMin);

    return cond;
}




#endif

