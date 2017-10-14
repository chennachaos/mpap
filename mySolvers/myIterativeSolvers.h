
#ifndef MY_ITERATIVE_SOLVERS_H
#define MY_ITERATIVE_SOLVERS_H


#include "headersBasic.h"
#include "myBLASroutines.h"
#include "omp.h"

using namespace std;
using namespace Eigen;



namespace myIterSolvers {




/*************************************************************

Conjugate gradient algorithm for the solution of the matrix system

  x = Ab

Input  --- Matrix A, rhs vector b
Output --- vector x

*************************************************************/


/** \internal Low-level bi conjugate gradient stabilized algorithm
  * \param mat The matrix A
  * \param rhs The right hand side vector b
  * \param x On input and initial solution, on output the computed solution.
  * \param precond A preconditioner being able to efficiently solve for an
  *                approximation of Ax=b (regardless of b)
  * \param iters On input the max number of iteration, on output the number of performed iterations.
  * \param tol_error On input the tolerance error, on output an estimation of the relative error.
  * \return false in the case of numerical issue, for example a break down of BiCGSTAB. 
  */


template<typename MatrixType, typename VectorType>
bool myPCG(const MatrixType& mat, const VectorType& rhs, VectorType& x, Index& Iters, typename VectorType::RealScalar& tol_error)
{
  using std::sqrt;
  using std::abs;

  typedef typename VectorType::RealScalar RealScalar;
  typedef typename VectorType::Scalar Scalar;

  RealScalar tol = tol_error;

  Index maxIters = Iters;

  Index size = mat.cols();

  VectorType p = VectorType::Zero(size);
  VectorType z = VectorType::Zero(size);
  VectorType v = VectorType::Zero(size);
  VectorType r = VectorType::Zero(size);

  //VectorType  diagVec = mat.diagonal().cwiseInverse();

  r = rhs - mat * x;
  z = r;
  //z = diagVec.cwiseProduct(r);
  p = z;
  
  RealScalar rhs_norm = rhs.norm();
  //RealScalar rhsNorm2 = rhs.squaredNorm();

  if(rhs_norm == 0.0)
  {
    x.setZero();
    return true;
  }

  Scalar rho    = 1.0;
  Scalar alpha  = 1.0;
  Scalar rho_old = 1.0;
  Scalar beta   = 1.0;
  RealScalar tmp;


  RealScalar tol2 = tol*rhs_norm;
  //RealScalar threshold = tol*tol*rhsNorm2;
  //RealScalar eps2 = NumTraits<Scalar>::epsilon()*NumTraits<Scalar>::epsilon();
  Index ii = 0;

  while( ii < maxIters )
  {
    rho = r.dot(z);

    v = mat*p;

    alpha = rho / r.dot(v);

    x += alpha * p;

    r -= alpha * v;

    tmp = r.norm();
    //tmp = r.squaredNorm();

    cout << ii << '\t' << tmp << endl;

    if( tmp < tol2 )
      break;

    //z = precond.solve(r);
    z = r;
    //z = diagVec.cwiseProduct(r);

    beta = r.dot(z)/rho;
    p = z + beta * p;

    ++ii;
  }
  r = rhs - mat * x;
  tol_error = r.norm()/rhs_norm;

  Iters = ii;

  return true;
}


template<typename MatrixType, typename VectorType>
bool myPCGserial(const MatrixType& mat, const VectorType& rhs, VectorType& x, Index& Iters, typename VectorType::RealScalar& tol_error)
{
  using std::sqrt;
  using std::abs;

  typedef typename VectorType::RealScalar RealScalar;
  typedef typename VectorType::Scalar Scalar;

  RealScalar tol = tol_error;

  Index maxIters = Iters;

//cout << mat << endl;

  Index size = mat.rows();

//cout << " AAAAAAAA " << endl;
//cout << size << endl;
//cout << " AAAAAAAA " << endl;

  VectorType p = VectorType::Zero(size);
  VectorType z = VectorType::Zero(size);
  VectorType v = VectorType::Zero(size);
  VectorType r = VectorType::Zero(size);
//cout << " CCCCCCCC " << endl;
  myMatVecMult(mat, x, p);

//cout << " AAAAAAAA " << endl;

  r = rhs - p;
  z = r;
  p = z;
  
  RealScalar rhs_norm = rhs.norm();

  if(rhs_norm == 0.0)
  {
    x.setZero();
    return true;
  }

  Scalar rho    = 1.0;
  Scalar alpha  = 1.0;
  Scalar rho_old = 1.0;
  Scalar beta   = 1.0;
  RealScalar tmp;


  RealScalar tol2 = tol*rhs_norm;
  //RealScalar eps2 = NumTraits<Scalar>::epsilon()*NumTraits<Scalar>::epsilon();
  Index ii = 0;

  while( ii < maxIters )
  {
    //rho = r.dot(z);
    rho = myDotProduct(&r(0), &z(0), size);

    //v = mat*p;
    myMatVecMult(mat, p, v);

    //alpha = rho / r.dot(v);
    alpha = rho / myDotProduct(&r(0), &v(0), size);

    //x += alpha * p;
    //x = alpha*p + x;
    myVecAXPBY(size, alpha, &p(0), 1.0, &x(0));

    //r -= alpha * v;
    //r = r - alpha * v;
    //r = - alpha * v + r;
    myVecAXPBY(size, -alpha, &v(0), 1.0, &r(0));

    //tmp = r.norm();
    tmp = sqrt(myDotProduct(&r(0), &r(0), size));

    cout << ii << '\t' << tmp << endl;

    if( tmp < tol2 )
      break;

    //z = precond.solve(r);
    z = r;

    //beta = r.dot(z)/rho;
    beta = myDotProduct(&r(0), &z(0), size)/rho;

    //p = z + beta * p;
    myVecAXPBY(size, 1.0, &z(0), beta, &p(0));

    ++ii;
  }

  tol_error = tmp/rhs_norm;

  Iters = ii;

  return true;
}


/** \internal Low-level bi conjugate gradient stabilized algorithm
  * \param mat The matrix A
  * \param rhs The right hand side vector b
  * \param x On input and initial solution, on output the computed solution.
  * \param precond A preconditioner being able to efficiently solve for an
  *                approximation of Ax=b (regardless of b)
  * \param iters On input the max number of iteration, on output the number of performed iterations.
  * \param tol_error On input the tolerance error, on output an estimation of the relative error.
  * \return false in the case of numerical issue, for example a break down of BiCGSTAB. 
  */


//template<typename MatrixType, typename Rhs, typename Dest, typename Preconditioner>
//bool bicgstab(const MatrixType& mat, const Rhs& rhs, Dest& x,
//              const Preconditioner& precond, Index& iters, typename Dest::RealScalar& tol_error)

template<typename MatrixType, typename VectorType>
bool myBiCGSTAB(const MatrixType& mat, const VectorType& rhs, VectorType& x, Index& Iters, typename VectorType::RealScalar& tol_error)
{
  using std::sqrt;
  using std::abs;

  typedef typename VectorType::RealScalar RealScalar;
  typedef typename VectorType::Scalar Scalar;

  RealScalar tol = tol_error;
  Index maxIters = Iters;

  Index n = mat.cols();
  VectorType r  = rhs - mat * x;
  VectorType r0 = r;
  
  RealScalar r0_sqnorm = r0.squaredNorm();
  RealScalar rhs_sqnorm = rhs.squaredNorm();
  if(rhs_sqnorm == 0)
  {
    x.setZero();
    return true;
  }
  Scalar rho    = 1.0;
  Scalar alpha  = 1.0;
  Scalar w      = 1.0;
  
  VectorType v = VectorType::Zero(n), p = VectorType::Zero(n);
  VectorType y(n),  z(n);
  VectorType kt(n), ks(n);

  VectorType s(n), t(n);

  RealScalar tol2 = tol*tol*rhs_sqnorm;
  RealScalar eps2 = NumTraits<Scalar>::epsilon()*NumTraits<Scalar>::epsilon();
  Index ii = 0;
  Index restarts = 0;

  while ( r.squaredNorm() > tol2 && ii<maxIters )
  {
    Scalar rho_old = rho;

    rho = r0.dot(r);
    if (abs(rho) < eps2*r0_sqnorm)
    {
      // The new residual vector became too orthogonal to the arbitrarily chosen direction r0
      // Let's restart with a new r0:
      r  = rhs - mat * x;
      r0 = r;
      rho = r0_sqnorm = r.squaredNorm();
      if(restarts++ == 0)
        ii = 0;
    }
    Scalar beta = (rho/rho_old) * (alpha / w);
    p = r + beta * (p - w * v);
    
    y = p;
    //y = precond.solve(p);
    
    v.noalias() = mat * y;

    alpha = rho / r0.dot(v);
    s = r - alpha * v;

    z = s;
    //z = precond.solve(s);
    t.noalias() = mat * z;

    RealScalar tmp = t.squaredNorm();
    if(tmp>RealScalar(0))
      w = t.dot(s) / tmp;
    else
      w = Scalar(0);
    x += alpha * y + w * z;
    r = s - w * t;

    ++ii;
  }

  tol_error = sqrt(r.squaredNorm()/rhs_sqnorm);

  Iters = ii;

  return true; 
}


} // end namespace myIterSolvers



#endif // MY_ITERATIVE_SOLVERS_H


