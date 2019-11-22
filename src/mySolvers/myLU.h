

#ifndef MY_LU_DECOMPOSITION_H
#define MY_LU_DECOMPOSITION_H


#include "headersBasic.h"

typedef  SparseMatrix<double, RowMajor>  SparseMatrixXd;


/**
** LU decomposition of matrix 'A'
** using Doolittle algorithm - 1s on the diagonal of L
** Lower part is stored in the matrix 'L'
** Upper part is stored in the matrix 'U'
** Input matrix 'A' is modified during the process
**/

void myLU_algo1(MatrixXd& A, MatrixXd& L, MatrixXd& U)
{
  // No row exchanges
  // Square matrix
  int k, i, j, n=A.rows();

  //# pragma omp parallel shared (a, n, L, U) private ( k, i, j )
  //{
    //# pragma omp for schedule ( static, 10 )
    for(i=0; i<n; i++)
    {
      L(i,i) = 1.0 ;
      for(k=i+1; k<n; k++)
      {
        L(k,i) = A(k,i)/A(i,i) ;

        for(j=i+1; j<n; j++)
        {
          A(k,j) = A(k,j) - L(k,i)*A(i,j);
        }
        //# pragma omp flush (A)
      }
      //# pragma omp nowait
      for(j=i; j<n; j++)
      {
        U(i,j) = A(i,j) ;
      }
    }
  //}
  return;
}




/**
** LU decomposition of matrix 'A'
** using Doolittle algorithm - 1s on the diagonal of L
** Lower part is stored in the matrix 'L'
** Upper part is stored in the matrix 'U'
** Input matrix 'A' does not get modified
**/

void myLU_algo2(MatrixXd& A, MatrixXd& L, MatrixXd& U)
{
  // No row exchanges
  // Square matrix
  int k, i, j, n=A.rows();

  for(i=0; i<n; i++)
  {
      L(i,i) = 1.0;

      //  if(j>i) L(i,j) = 0.0

      for(j=0; j<i; j++)
      {
        L(i,j) = A(i,j);

        for(k=0; k<j; k++)
        {
          L(i,j) = L(i,j) - U(k,j)*L(i,k) ;
        }

        L(i,j) = L(i,j) / U(j,j);
      }

      //if(j<i) U(i,j) = 0.0;

      for(j=i; j<n; j++)
      {
        U(i,j) = A(i,j);

        for(k=0; k<i; k++)
        {
          U(i,j) = U(i,j) - U(k,j)*L(i,k) ;
        }
      }
    }

  return;
}


/**
** LU decomposition of matrix 'A'
** using Doolittle algorithm - 1s on the diagonal of L
** L and U are stored in the matrix 'mLU'
** For solving the lower part 'L' needs to be appendid with 1s on the diagonal
**/

void myLU_algo3(MatrixXd& A, MatrixXd& mLU)
{
  // No row exchanges
  // Square matrix
  int k, i, j, n=A.rows();

    for(i=0; i<n; i++)
    {
        //L(i,i) = 1.0;

        //  if(j>i) L(i,j) = 0.0

        for(j=0; j<i; j++)
        {
           mLU(i,j) = A(i,j);
           for(k=0; k<j; k++)
           {
              mLU(i,j) = mLU(i,j) - mLU(k,j)*mLU(i,k) ;
           }
           mLU(i,j) = mLU(i,j) / mLU(j,j);
        }

        //if(j<i) U(i,j) = 0.0;

        for(j=i; j<n; j++)
        {
           mLU(i,j) = A(i,j);
           for(k=0; k<i; k++)
           {
              mLU(i,j) = mLU(i,j) - mLU(k,j)*mLU(i,k) ;
           }
        }

    }

  return;
}





/**
** LU decomposition of matrix 'A'

** using Crout's algorithm - 1s on the diagonal of U

** Lower part is stored in the matrix 'L'
** Upper part is stored in the matrix 'U'
** Input matrix 'A' does not get modified
**/

void myLU_algo4(MatrixXd& A, MatrixXd& L, MatrixXd& U)
{
  // No row exchanges
  // Square matrix
  int k, i, j, n=A.rows();

    for(i=0; i<n; i++)
    {
      for(j=i; j<n; j++)
      {
        L(j,i) = A(j,i);
        for(k=0; k<i; k++)
        {
          L(j,i) = L(j,i) - L(j,k)*U(k,i);
        }
      }

      U(i,i) = 1.0;
      for(j=i+1; j<n; j++)
      {
        U(i,j) = A(i,j) / L(i,i);
        for(k=0; k<i; k++)
        {
          U(i,j) = U(i,j) - ((L(i,k)*U(k,j)) / L(i,i));
        }
      }
    }

    return;
}



/**
** LU decomposition of matrix 'A'
** using Doolittle algorithm - 1s on the diagonal of L
** Lower part is stored in the matrix 'L'
** Upper part is stored in the matrix 'U'
** Input matrix 'A' does not get modified
**/

void myLU_algo2(SparseMatrixXd& A, SparseMatrixXd& L, SparseMatrixXd& U)
{
  // No row exchanges
  // Square matrix
  int k, i, j, n=A.rows();
  double  val=0.0, temp=0.0;
  int k1, ik, ikk, k2;

  const int* rowPtrA   = A.outerIndexPtr();
  const int* colIndxA  = A.innerIndexPtr();
  const double* arrayA = A.valuePtr();

  i = A.nonZeros()*20;

  L.resize(n,n);
  L.reserve(i);

  U.resize(n,n);
  U.reserve(i);

  VectorXd u(n) ;   // real values of the row -- maximum size is n --
  VectorXi ju(n);   // column position of the values in u -- maximum size  is n
  VectorXi jr(n);   // Indicate the position of the nonzero elements in the vector u -- A zero location is indicated by -1

  // Initialization
  jr.fill(-1);
  ju.fill(0);
  u.fill(0);


  for(i=0; i<n; i++)
  {
      L.coeffRef(i, i) = 1.0;

      //  if(j>i) L(i,j) = 0.0

      k1  = rowPtrA[i];
      ikk = (rowPtrA[i+1]-rowPtrA[i]);
      //cout << ii << '\t' << k1 << '\t' << ikk << endl;       printf("\n\n");

      /*
       // Iterate through the current row ii
      for(ik=0; ik<(ikk-1); ik++)
      {
        k2 = k1+ik;
        j = colIndxA[k2];
        //cout << " inside " << k2 << '\t' << k << '\t' << temp << endl;

        if(j<i)
        {
          val = arrayA[k2];

          for(k=0; k<j; k++)
          {
            val = val - U(k,j)*L(i,k) ;
          }

          L.coeffRef(i,j) = val / u(j);
        }
      }
      */


      //if(j<i) U(i,j) = 0.0;

      /*
      for(j=i; j<n; j++)
      {
        U(i,j) = A(i,j);

        for(k=0; k<i; k++)
        {
          U(i,j) = U(i,j) - U(k,j)*L(i,k) ;
        }
      }
      */
    }

  L.finalize();
  L.makeCompressed();

  U.finalize();
  U.makeCompressed();

  return;
}




//!  This function is used to find the solution to a system of equations,
/*!   A x = b, after LU decomposition of A has been found.
*    Within this routine, the elements of b are rearranged in the same way
*    that the rows of a were interchanged, using the order vector.
*    The solution is returned in x.
*
*
*  \param  a     - the LU decomposition of the original coefficient Matrix.
*  \param  b     - the vector of right-hand sides
*  \param       x     - the solution vector
*  \param    order - integer array of row order as arranged during pivoting
*
*/
void solvLU(const MatrixXd& L, const MatrixXd& U, const VectorXd& b, const VectorXi& order, VectorXd& x)
{
  int i,j,n = L.rows();
  double sum;

    /* rearrange the elements of the b vector. x is used to hold them. */

    for(i=0; i<n; i++)
    {
      x[i] = b[order[i]];
    }

    // Forward solve Ly = b

    for(i=0; i<n; i++)
    {
      for (j=0; j<i; j++)
        x[i] -= L(i,j)*x[j];

      x[i] /= L(i,i);
    }

    // Backward solve Ux = y

    for(i=n-1; i>=0; i--)
    {
      for(j=i+1; j<n; j++) 
        x[i] -= U(i,j) * x[j];

      x[i] /= U(i,i);
    }

    return;
}


void solvLU(const SparseMatrixXd& L, const SparseMatrixXd& U, const VectorXd& b, const VectorXi& order, VectorXd& x)
{
  int ii,j, k, n = L.rows();
  double  temp;
  int k1, ik, ikk, k2;

    /* rearrange the elements of the b vector. x is used to hold them. */

    for(ii=0; ii<n; ii++)
    {
      x[ii] = b[ii];
      //x[i] = b[order[i]];
    }

    // Forward solve Ly = b

    int* rowPtrL   = L.outerIndexPtr();
    int* colIndxL  = L.innerIndexPtr();
    double* arrayL = L.valuePtr();

    for(ii=0; ii<n; ii++)
    {
      k1  = rowPtrL[ii];
      ikk = (rowPtrL[ii+1]-rowPtrL[ii]);
      //cout << ii << '\t' << k1 << '\t' << ikk << endl;       printf("\n\n");

       // Iterate through the current row ii
      //for(SparseMatrixXd::InnerIterator it(L, ii); it; ++it)
      for(ik=0; ik<(ikk-1); ik++)
      {
        k2 = k1+ik;
        k = colIndxL[k2];
        temp = arrayL[k2];
        //cout << " inside " << k2 << '\t' << k << '\t' << temp << endl;

        x[ii] -= temp*x[k];
      }

      x[ii] /=  arrayL[k1+ik];
    }

    //printVector(x);    printf("\n\n");

    // Backward solve Ux = y

    int* rowPtrU   = U.outerIndexPtr();
    int* colIndxU  = U.innerIndexPtr();
    double* arrayU = U.valuePtr();

    for(ii=n-1; ii>=0; ii--)
    {
      k1  = rowPtrU[ii];
      ikk = (rowPtrU[ii+1]-rowPtrU[ii]);

      //cout << ii << '\t' << k1 << '\t' << ikk << endl;       printf("\n\n");

       // Iterate through the current row ii
      //for(SparseMatrixXd::InnerIterator it(U, ii); it; ++it)
      for(ik=(ikk-1); ik>0; ik--)
      {
        k2 = k1+ik;
        k = colIndxU[k2];
        temp = arrayU[k2];

        //cout << " inside " << k2 << '\t' << k << '\t' << temp << endl;

        x[ii] -= temp*x[k];
      }

      x[ii] /=  arrayU[k1];
    }

    return;
}



#endif







