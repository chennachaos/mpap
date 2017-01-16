
#ifndef MY_BLAS_ROUTINES_H
#define MY_BLAS_ROUTINES_H


#include "headersBasic.h"
#include "omp.h"

using namespace std;
using namespace Eigen;



namespace myIterSolvers {


//template<typename myMat, typename myVec>
//void myMatVecMult(myMat& A, myVec& x, myVec& y)

//void myMatVecMult(const int nRow, const int* csr, const int* col, const double* val, const double* x, double* y)

//template<typename ScalarType>
//double myDotProduct(ScalarType* x, ScalarType* y, int size)

//template<typename ScalarType>
//void myVecAdd2(int size, ScalarType* x, ScalarType* y, ScalarType* w)

//template<typename ScalarType>
//void myVecAdd3(int size, ScalarType* x, ScalarType* y, ScalarType* z, ScalarType* w)

//template<typename ScalarType>
//void myVecAXPBY(int size, ScalarType a, ScalarType* x, ScalarType b, ScalarType* y)

//template<typename ScalarType>
//void myVecZAXPBY(int size, ScalarType a, ScalarType* x, ScalarType b, ScalarType* y, ScalarType* z)


typedef int Index;


/*************************************************************

subroutine to multiply a matrix (A) with a vector (x)

  y = Ax

Input  --- matrix A, vector x
Output --- vector y

*************************************************************/

template<typename myMat, typename myVec>
void myMatVecMult(myMat& A, myVec& x, myVec& y)
{
    int ii, jj, size = x.rows();
    double  val;

    const int *csr  = A.outerIndexPtr();
    const int *col  = A.innerIndexPtr();
    const double *data = A.valuePtr();

    //#pragma omp parallel for  shared(size, csr, col, val) private(ii, jj) 
    #pragma omp parallel for  shared(size, csr, col, val) private(ii, jj) 
    for(ii=0; ii<size; ii++)
    {
       val = 0.0;
       for(jj=csr[ii]; jj<csr[ii+1]; jj++)
         val += data[jj] * x[col[jj]];

       y[ii] = val;
    }

    return;
}


/*************************************************************

subroutine to multiply a matrix (A) with a vector (x)
Matrix A is assumed to be stored in CSR format 
and corresponding details are supplied as input

  y = Ax

Input  --- csr of A, column pointer of A, data array of A, vector x
           rows in A, cols in A
Output --- vector y

*************************************************************/

void myMatVecMult(const int nRow, const int* csr, const int* col, const double* data, const double* x, double* y)
{
    int ii, jj;
    double  val;

    for(ii=0; ii<nRow; ii++)
    {
       val = 0.0;
       //#pragma omp parallel for  shared(size, x, y) private(ii)
       for(jj=csr[ii]; jj<csr[ii+1]; jj++)
         val += data[jj] * x[col[jj]];

       y[ii] = val;
    }

    return;
}

/*************************************************************

subroutine to compute dot product of two vectors (x and y) 

  val = x.y = x^T y = y^T x

Input  --- vector x, vector y
Output --- x.y

*************************************************************/


template<typename ScalarType>
double myDotProduct(ScalarType* x, ScalarType* y, int size)
{
    //double  val=0.0;
    //for(int ii=0; ii<size; ii++)
      //val += x[ii] * y[ii];

    int ii;
    double  val=0.0;

    #pragma omp parallel for  shared(size, x, y) private(ii) reduction(+:val)
    for(ii=0; ii<size; ii++)
    {
       val += x[ii] * y[ii];
    }

    return val;
}


/*************************************************************

subroutine to compute

  w = x + y

Input  --- vector x, vector y
Output --- vector y

*************************************************************/

template<typename ScalarType>
void myVecAdd2(int size, ScalarType* x, ScalarType* y, ScalarType* w)
{
    int ii;

    #pragma omp parallel for shared(size, x, y, w) private(ii)
    for(ii=0; ii<size; ii++)
    {
       w[ii] = x[ii] + y[ii];
    }

    return;
}


/*************************************************************

subroutine to compute

  w = x + y + z

Input  --- vector x, vector y, vector z
Output --- vector w

*************************************************************/

template<typename ScalarType>
void myVecAdd3(int size, ScalarType* x, ScalarType* y, ScalarType* z, ScalarType* w)
{
    int ii;

    #pragma omp parallel for shared(size, x, y, z, w) private(ii)
    for(ii=0; ii<size; ii++)
    {
       w[ii] = x[ii] + y[ii] + z[ii];
    }

    return;
}

/*************************************************************

subroutine to scale a vector by a scalar (alpha)

  x = alpha*x

Input  --- vector x, scalar alpha
Output --- vector x

*************************************************************/


template<typename ScalarType>
void myVecScale(int size, ScalarType* x, ScalarType alpha)
{
    int ii;

    #pragma omp parallel for shared(size, x, alpha) private(ii)
    for(ii=0; ii<size; ii++)
    {
       x[ii] = x[ii] * alpha;
    }

    return;
}


/*************************************************************

subroutine to compute

  y = a*x + b*y

Input  --- vector x, vector y, scalar a, scalar b
Output --- vector y

*************************************************************/


template<typename ScalarType>
void myVecAXPBY(int size, ScalarType a, ScalarType* x, ScalarType b, ScalarType* y)
{
    int ii;

    #pragma omp parallel for shared(size, a, x, b, y) private(ii)
    for(ii=0; ii<size; ii++)
    {
       y[ii] = a*x[ii] + b*y[ii];
    }

    return;
}

/*************************************************************

subroutine to compute

  z = a*x + b*y

Input  --- vector x, vector y, scalar a, scalar b
Output --- vector z

*************************************************************/


template<typename ScalarType>
void myVecZAXPBY(int size, ScalarType a, ScalarType* x, ScalarType b, ScalarType* y, ScalarType* z)
{
    int ii;
    #pragma omp parallel for shared(size, a, x, b, y, z) private(ii)
    for(ii=0; ii<size; ii++)
    {
       z[ii] = a*x[ii] + b*y[ii];
    }

    return;
}



} // end namespace myIterSolvers



#endif // MY_BLAS_ROUTINES_H


