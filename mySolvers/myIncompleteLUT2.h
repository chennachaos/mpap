

#ifndef EIGEN_MY_INCOMPLETE_LUT2_H
#define EIGEN_MY_INCOMPLETE_LUT2_H

#include <vector>

#include "headersEigen.h"

#include "Eigen/src/IterativeLinearSolvers/IncompleteLUT.h"

#include "Solver.h"
#include "FunctionsSolver.h"

using std::vector;

using namespace Eigen;
using namespace Eigen::internal;


namespace myIterSolvers
{


template <typename _Scalar, typename _StorageIndex = int>
class  myIncompleteLUT2 : public SparseSolverBase<myIncompleteLUT2<_Scalar, _StorageIndex> >
{
  protected:
    typedef SparseSolverBase<myIncompleteLUT2> Base;
    using Base::m_isInitialized;
  public:
    typedef _Scalar Scalar;
    typedef _StorageIndex StorageIndex;
    typedef typename NumTraits<Scalar>::Real RealScalar;
    typedef Matrix<Scalar,Dynamic,1> Vector;
    typedef Matrix<StorageIndex,Dynamic,1> VectorI;
    //typedef SparseMatrix<Scalar,RowMajor,StorageIndex> FactorType;
    typedef SparseMatrix<Scalar,RowMajor> FactorType;

    void   *PT[64];

    Scalar  DPARM[64], *array;

    vector<StorageIndex>  csr, col;
    vector<StorageIndex>  perm;

    StorageIndex    SOLVER, MTYPE, MAXFCT, MNUM, MSGLVL, NRHS, IPARM[64];
    StorageIndex   nRow, nCol;
    int  currentStatus;
    VectorXd  soln;

  public:
    
    // this typedef is only to export the scalar type and compile-time dimensions to solve_retval
    typedef Matrix<Scalar,Dynamic,Dynamic> MatrixType;
    
    myIncompleteLUT2()
      : m_droptol(NumTraits<Scalar>::dummy_precision()), m_fillfactor(10),
        m_analysisIsOk(false), m_factorizationIsOk(false), nRow(0)
    {}

    ~myIncompleteLUT2();

    int initialise();
   
    int factorise();

    VectorXd&  solve(const VectorXd&  rhsVec);

    VectorXd&  factoriseAndSolve(const VectorXd&  rhsVec);

    void free();


    template<typename MatrixType>
    explicit  myIncompleteLUT2(const MatrixType& mat, const RealScalar& droptol=NumTraits<Scalar>::dummy_precision(), int fillfactor = 10)
      : m_droptol(droptol),m_fillfactor(fillfactor),
        m_analysisIsOk(false),m_factorizationIsOk(false)
    {
      eigen_assert(fillfactor != 0);
      compute(mat); 
    }
    
    Index rows() const { return m_lu.rows(); }
    
    Index cols() const { return m_lu.cols(); }

    /** \brief Reports whether previous computation was successful.
      *
      * \returns \c Success if computation was succesful,
      *          \c NumericalIssue if the matrix.appears to be negative.
      */
    ComputationInfo info() const
    {
      eigen_assert(m_isInitialized && "myIncompleteLUT is not initialized.");
      return m_info;
    }
    
    template<typename MatrixType>
    void analyzePattern(const MatrixType& amat);
    
    template<typename MatrixType>
    void factorize(const MatrixType& amat);
    
    /**
      * Compute an incomplete LU factorization with dual threshold on the matrix mat
      * No pivoting is done in this version
      * 
      **/
    template<typename MatrixType>
    myIncompleteLUT2&  compute(const MatrixType& amat)
    {
      m_analysisIsOk = true;
      m_factorizationIsOk = false;
      m_isInitialized = false;

      //analyzePattern(amat);
      //factorise();

      factorize(amat);

      initialise();
      
      return *this;
    }

    void setDroptol(const RealScalar& droptol); 
    void setFillfactor(int fillfactor); 

    template<typename Rhs, typename Dest>
    void _solve_impl(const Rhs& b, Dest& x) const
    {
      x = m_Pinv * b;  
      x = m_lu.template triangularView<UnitLower>().solve(x);
      x = m_lu.template triangularView<Upper>().solve(x);
      x = m_P * x; 
    }



protected:

    /** keeps off-diagonal entries; drops diagonal entries */
    struct keep_diag {
      inline bool operator() (const Index& row, const Index& col, const Scalar&) const
      {
        return row!=col;
      }
    };


protected:

    FactorType m_lu;
    RealScalar m_droptol;
    int m_fillfactor;
    bool m_analysisIsOk;
    bool m_factorizationIsOk;
    ComputationInfo m_info;
    PermutationMatrix<Dynamic,Dynamic,StorageIndex> m_P;     // Fill-reducing permutation
    PermutationMatrix<Dynamic,Dynamic,StorageIndex> m_Pinv;  // Inverse permutation

};



/**
 * destructor 
 **/ 
template<typename Scalar, typename StorageIndex>
myIncompleteLUT2<Scalar,StorageIndex>::~myIncompleteLUT2()
{
  int phase = -1, error = 0;
  double ddum;

  //pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
    //       &nRow, array, &csr[0], &col[0], &perm[0], &NRHS,
      //     IPARM, &MSGLVL, &ddum, &ddum, &error, DPARM);
  if(nRow != 0)
    free();

}


/**
 * Set control parameter droptol
 *  \param droptol   Drop any element whose magnitude is less than this tolerance 
 **/ 
template<typename Scalar, typename StorageIndex>
void myIncompleteLUT2<Scalar,StorageIndex>::setDroptol(const RealScalar& droptol)
{
  this->m_droptol = droptol;   
}



/**
 * Set control parameter fillfactor
 * \param fillfactor  This is used to compute the  number @p fill_in of largest elements to keep on each row. 
 **/ 
template<typename Scalar, typename StorageIndex>
void myIncompleteLUT2<Scalar,StorageIndex>::setFillfactor(int fillfactor)
{
  this->m_fillfactor = fillfactor;   
}



template <typename Scalar, typename StorageIndex>
template<typename _MatrixType>
void myIncompleteLUT2<Scalar,StorageIndex>::analyzePattern(const _MatrixType& amat)
{
  // Compute the Fill-reducing permutation
  SparseMatrix<Scalar,ColMajor, StorageIndex> mat1 = amat;
  SparseMatrix<Scalar,ColMajor, StorageIndex> mat2 = amat.transpose();
  // Symmetrize the pattern
  // FIXME for a matrix with nearly symmetric pattern, mat2+mat1 is the appropriate choice.
  //       on the other hand for a really non-symmetric pattern, mat2*mat1 should be prefered...
  SparseMatrix<Scalar,ColMajor, StorageIndex> AtA = mat2 + mat1;
  AtA.prune(keep_diag());
  internal::minimum_degree_ordering<Scalar, StorageIndex>(AtA, m_P);  // Then compute the AMD ordering...

  m_Pinv  = m_P.inverse(); // ... and the inverse permutation

  m_analysisIsOk = true;
  m_factorizationIsOk = false;
  m_isInitialized = false;
}



template <typename Scalar, typename StorageIndex>
template<typename _MatrixType>
void  myIncompleteLUT2<Scalar,StorageIndex>::factorize(const _MatrixType& amat)
{
/*
  using std::sqrt;
  using std::swap;
  using std::abs;
  using internal::convert_index;

  eigen_assert((amat.rows() == amat.cols()) && "The factorization should be done on a square matrix");
  Index n = amat.cols();  // Size of the matrix
  m_lu.resize(n,n);
  // Declare Working vectors and variables
  Vector u(n) ;     // real values of the row -- maximum size is n --
  VectorI ju(n);   // column position of the values in u -- maximum size  is n
  VectorI jr(n);   // Indicate the position of the nonzero elements in the vector u -- A zero location is indicated by -1

  // Apply the fill-reducing permutation
  eigen_assert(m_analysisIsOk && "You must first call analyzePattern()");
  SparseMatrix<Scalar,RowMajor, StorageIndex> mat;
  mat = amat.twistedBy(m_Pinv);

  // Initialization
  jr.fill(-1);
  ju.fill(0);
  u.fill(0);

  // number of largest elements to keep in each row:
  Index fill_in = (amat.nonZeros()*m_fillfactor)/n + 1;
  if (fill_in > n) fill_in = n;

  // number of largest nonzero elements to keep in the L and the U part of the current row:
  Index nnzL = fill_in/2;
  Index nnzU = nnzL;
  m_lu.reserve(n * (nnzL + nnzU + 1));

  // global loop over the rows of the sparse matrix
  for (Index ii = 0; ii < n; ii++)
  {
    // 1 - copy the lower and the upper part of the row i of mat in the working vector u

    Index sizeu = 1; // number of nonzero elements in the upper part of the current row
    Index sizel = 0; // number of nonzero elements in the lower part of the current row
    ju(ii)    = convert_index<StorageIndex>(ii);
    u(ii)     = 0;
    jr(ii)    = convert_index<StorageIndex>(ii);
    RealScalar rownorm = 0;

    typename FactorType::InnerIterator j_it(mat, ii); // Iterate through the current row ii
    for (; j_it; ++j_it)
    {
      Index k = j_it.index();
      if (k < ii)
      {
        // copy the lower part
        ju(sizel) = convert_index<StorageIndex>(k);
        u(sizel) = j_it.value();
        jr(k) = convert_index<StorageIndex>(sizel);
        ++sizel;
      }
      else if (k == ii)
      {
        u(ii) = j_it.value();
      }
      else
      {
        // copy the upper part
        Index jpos = ii + sizeu;
        ju(jpos) = convert_index<StorageIndex>(k);
        u(jpos) = j_it.value();
        jr(k) = convert_index<StorageIndex>(jpos);
        ++sizeu;
      }
      rownorm += numext::abs2(j_it.value());
    }

    // 2 - detect possible zero row
    if(rownorm==0)
    {
      m_info = NumericalIssue;
      return;
    }
    // Take the 2-norm of the current row as a relative tolerance
    rownorm = sqrt(rownorm);

    // 3 - eliminate the previous nonzero rows
    Index jj = 0;
    Index len = 0;
    while (jj < sizel)
    {
      // In order to eliminate in the correct order,
      // we must select first the smallest column index among  ju(jj:sizel)
      Index k;
      Index minrow = ju.segment(jj,sizel-jj).minCoeff(&k); // k is relative to the segment
      k += jj;
      if (minrow != ju(jj))
      {
        // swap the two locations
        Index j = ju(jj);
        swap(ju(jj), ju(k));
        jr(minrow) = convert_index<StorageIndex>(jj);
        jr(j) = convert_index<StorageIndex>(k);
        swap(u(jj), u(k));
      }
      // Reset this location
      jr(minrow) = -1;

      // Start elimination
      typename FactorType::InnerIterator ki_it(m_lu, minrow);
      while (ki_it && ki_it.index() < minrow) ++ki_it;
      eigen_internal_assert(ki_it && ki_it.col()==minrow);
      Scalar fact = u(jj) / ki_it.value();

      // drop too small elements
      if(abs(fact) <= m_droptol)
      {
        jj++;
        continue;
      }

      // linear combination of the current row ii and the row minrow
      ++ki_it;
      for (; ki_it; ++ki_it)
      {
        Scalar prod = fact * ki_it.value();
        Index j     = ki_it.index();
        Index jpos  = jr(j);
        if (jpos == -1) // fill-in element
        {
          Index newpos;
          if (j >= ii) // dealing with the upper part
          {
            newpos = ii + sizeu;
            sizeu++;
            eigen_internal_assert(sizeu<=n);
          }
          else // dealing with the lower part
          {
            newpos = sizel;
            sizel++;
            eigen_internal_assert(sizel<=ii);
          }
          ju(newpos) = convert_index<StorageIndex>(j);
          u(newpos) = -prod;
          jr(j) = convert_index<StorageIndex>(newpos);
        }
        else
          u(jpos) -= prod;
      }
      // store the pivot element
      u(len)  = fact;
      ju(len) = convert_index<StorageIndex>(minrow);
      ++len;

      jj++;
    } // end of the elimination on the row ii

    // reset the upper part of the pointer jr to zero
    for(Index k = 0; k <sizeu; k++) jr(ju(ii+k)) = -1;

    // 4 - partially sort and insert the elements in the m_lu matrix

    // sort the L-part of the row
    sizel = len;
    len = (std::min)(sizel, nnzL);
    typename Vector::SegmentReturnType ul(u.segment(0, sizel));
    typename VectorI::SegmentReturnType jul(ju.segment(0, sizel));
    internal::QuickSplit(ul, jul, len);

    // store the largest m_fill elements of the L part
    m_lu.startVec(ii);
    for(Index k = 0; k < len; k++)
      m_lu.insertBackByOuterInnerUnordered(ii,ju(k)) = u(k);

    // store the diagonal element
    // apply a shifting rule to avoid zero pivots (we are doing an incomplete factorization)
    if (u(ii) == Scalar(0))
      u(ii) = sqrt(m_droptol) * rownorm;
    m_lu.insertBackByOuterInnerUnordered(ii, ii) = u(ii);

    // sort the U-part of the row
    // apply the dropping rule first
    len = 0;
    for(Index k = 1; k < sizeu; k++)
    {
      if(abs(u(ii+k)) > m_droptol * rownorm )
      {
        ++len;
        u(ii + len)  = u(ii + k);
        ju(ii + len) = ju(ii + k);
      }
    }
    sizeu = len + 1; // +1 to take into account the diagonal element
    len = (std::min)(sizeu, nnzU);
    typename Vector::SegmentReturnType uu(u.segment(ii+1, sizeu-1));
    typename VectorI::SegmentReturnType juu(ju.segment(ii+1, sizeu-1));
    internal::QuickSplit(uu, juu, len);

    // store the largest elements of the U part
    for(Index k = ii + 1; k < ii + len; k++)
      m_lu.insertBackByOuterInnerUnordered(ii,ju(k)) = u(k);
  }
*/

  //
  Index n = amat.cols();  // Size of the matrix
  nRow = nCol = n;
  m_lu.resize(n,n);

  // number of largest nonzero elements to keep in the L and the U part of the current row:
  //Index  nnzLU = m_fillfactor*n;
  Index  nnzLU = amat.nonZeros()/4;
  Index  ii, jj, kk, ll, k1, k2;

  m_lu.reserve(nnzLU);
  
  Scalar  fact;

  // global loop over the rows of the sparse matrix
  for(ii=0; ii<n; ii++)
  {
    k1 = ii-m_fillfactor;
    k2 = ii+m_fillfactor;

    for(typename FactorType::InnerIterator j_it(amat, ii); j_it; ++j_it)
    {
      if(j_it && (j_it.col() >= k1) && (j_it.col() <= k2) )
      {
        //if( j_it.value()!=Scalar(0) )
          //fact = j_it.value();
        //else
          //fact = 1.0;

        m_lu.coeffRef(ii, j_it.col()) = j_it.value();
      }
    } // end of the elimination on the row ii
  }
  //

  
  //////////////////////////////////////////////////
  /*

  Index  nnzLU = amat.nonZeros()/4;
  Index  n = amat.rows();
  

  m_lu.resize(n,n);
  m_lu.reserve(nnzLU);

  RealScalar rownorm = 0.0;
  Index  k, ii;

  // global loop over the rows of the sparse matrix
  for(ii=0; ii<n; ii++)
  {
    rownorm = 0.0;

    //typename FactorType::InnerIterator j_it(mat, ii); // Iterate through the current row ii

    for(typename FactorType::InnerIterator j_it(amat, ii); j_it; ++j_it)
    {
      rownorm += numext::abs2(j_it.value());
    }

    // 2 - detect possible zero row
    if( abs(rownorm) <= 1.0e-12)
    {
      m_info = NumericalIssue;
      cerr << " myIncompleteLUT .... Numerical issues ... rownorm=0 " << endl;
      return;
    }
    // Take the 2-norm of the current row as a relative tolerance
    rownorm = sqrt(rownorm);

    // drop too small elements
    for(typename FactorType::InnerIterator j_it(amat, ii); j_it; ++j_it)
    {
      //rownorm += numext::abs2(j_it.value());

      if(abs(j_it.value()) >= m_droptol*rownorm)
      {
        m_lu.coeffRef(ii, j_it.col()) = j_it.value();
      }
    } // end of the elimination on the row ii
  }
  */
  cout << " AAAAAAAAAAAAAA " << endl;

  m_lu.finalize();
  m_lu.makeCompressed();

  m_factorizationIsOk = true;
  m_isInitialized = m_factorizationIsOk;
  m_info = Success;
}





template<typename Scalar, typename StorageIndex>
int  myIncompleteLUT2<Scalar,StorageIndex>::initialise()
{
  char fct[] = "myIncompleteLUT::initialise";

  nRow = nCol = m_lu.rows();

  cout << " nRow " << nRow << endl;

  soln.resize(nRow);
  soln.setZero();

  //if (currentStatus != PATTERN_OK)
    //{ prgWarning(1,fct,"prepare matrix pattern first!"); return 1; }

  double ddum;

  int  phase = 12, error = 0, mtxType = 11, idum;

  char *tmp;

  mtxType =  1;  // real and structurally symmetric
  mtxType = 11;  // real and unsymmetric

  SOLVER = 0;       // sparse direct solver
  MTYPE  = mtxType; // matrix type
  MAXFCT = 1;       // maximum number of factorisations of same sparsity pattern

  MNUM   = 1;       // which factorisation to use
  NRHS   = 1;       // number of right hand sides
  MSGLVL = 0;       // output message level (1 -> print statistical information)

  IPARM[0] = 0;     // PARADISO will set IPARM to default values

  //cout << " numProc " << numProc << endl;
  
  //IPARM[9] = 1;

  tmp = getenv("OMP_NUM_THREADS");

  if(tmp != NULL)
  {
    sscanf(tmp,"%d", &idum);
    //if (idum != IPARM[2]) prgError(1,fct,"set environment variable OMP_NUM_THREADS to numProc!");
  }
  else
  {
    idum = 1;
    //prgError(2,fct,"set environment variable OMP_NUM_THREADS!");
  }

  IPARM[2] = idum;  // number of processors (no default value available)


  pardisoinit_(PT, &MTYPE, &SOLVER, IPARM, DPARM, &error);

  if (error != 0)
  {
    if (error == -10) prgError(1,fct,"no license file found.");
    if (error == -11) prgError(2,fct,"license is expired.");
    if (error == -12) prgError(3,fct,"wrong username or hostname.");
  }

  cout << "\n\n PARDISO license check was successful  \n " << endl;

  int  *c1, *c2, ii;
  
  csr.resize(nRow+1);
  col.resize(m_lu.nonZeros());

  c1 = m_lu.outerIndexPtr();
  c2 = m_lu.innerIndexPtr();

  for(ii=0;ii<=nRow;ii++)
    csr[ii] = c1[ii] + 1;

  for(ii=0;ii<m_lu.nonZeros();ii++)
    col[ii] = c2[ii] + 1;

  array = m_lu.valuePtr();

  //cout << " hhhhhhhhhhhhhhhhh " << endl;
  cout << " ILU nnz = " << m_lu.nonZeros() << endl;

  //IPARM[4] = 0;  // user input permutation
  //IPARM[4] = 2;  // return the permutation

  perm.resize(nRow);

  //cout << " ppppppppppp " << IPARM[14] << '\t' << IPARM[15] <<  endl;

  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &nRow, array, &csr[0], &col[0], &perm[0], &NRHS,
           IPARM, &MSGLVL, &ddum, &soln[0], &error, DPARM);

  //cout << " llllllllll " << endl;

  if (error != 0)
  {
    COUT << "myIncompleteLUT  = " << error << "\n\n";
    prgError(4,fct,"symbolic factorisation failed.");
  }

  IPARM[5] = 0; // do not overwrite RHS with solution
  IPARM[7] = 1; // max number of iterative refinement steps

  currentStatus = INIT_OK;
  currentStatus = ASSEMBLY_OK;
  currentStatus = FACTORISE_OK;
  
  return 1;
}



template<typename Scalar, typename StorageIndex>
int  myIncompleteLUT2<Scalar,StorageIndex>::factorise()
{
  char fct[] = "myIncompleteLUT::factorise";
  
  //cout << " myIncompleteLUT::factorise " << endl;

  if (currentStatus != ASSEMBLY_OK) 
    { prgWarning(1,fct,"assemble matrix first!"); return 1; }

  int phase = 22, idum, error = 0;
  double ddum;

  //computerTime.go(fct);

  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &nRow, array, &csr[0], &col[0], &perm[0], &NRHS,
           IPARM, &MSGLVL, &ddum, &ddum, &error, DPARM);

  //solverTime.total     -= solverTime.factorise;
  //solverTime.factorise += computerTime.stop(fct);
  //solverTime.total     += solverTime.factorise;

  currentStatus = FACTORISE_OK;
  
  return 0;
}





template<typename Scalar, typename StorageIndex>
VectorXd&  myIncompleteLUT2<Scalar,StorageIndex>::solve(const VectorXd&  rhsVec)
{ 
  char fct[] = "myIncompleteLUT::solve";

  //cout << " myIncompleteLUT::solve " << endl;

  if (currentStatus != FACTORISE_OK)
  {
    prgWarning(1,fct,"factorise matrix first!");
  }

  int phase = 33, error = 0;
  soln.setZero();

  //computerTime.go(fct);

  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &nRow, array, &csr[0], &col[0], &perm[0], &NRHS,
           IPARM, &MSGLVL, &rhsVec[0], &soln[0], &error, DPARM);

  //solverTime.total -= solverTime.solve;
  //solverTime.solve += computerTime.stop(fct);
  //solverTime.total += solverTime.solve;

  return soln;
}



template<typename Scalar, typename StorageIndex>
VectorXd&  myIncompleteLUT2<Scalar,StorageIndex>::factoriseAndSolve(const VectorXd&  rhsVec)
{
  char fct[] = "myIncompleteLUT::factoriseAndSolve";

  if (currentStatus != ASSEMBLY_OK) 
    { prgWarning(1,fct,"assemble matrix first!"); return 1; }

  int phase = 23, error = 0;
  soln.setZero();
  //computerTime.go(fct);

  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &nRow, array, &csr[0], &col[0], &perm[0], &NRHS,
           IPARM, &MSGLVL, &rhsVec[0], &soln[0], &error, DPARM);
  
  //solverTime.total             -= solverTime.factoriseAndSolve;
  //solverTime.factoriseAndSolve += computerTime.stop(fct);
  //solverTime.total             += solverTime.factoriseAndSolve;

  currentStatus = FACTORISE_OK;
  
  return soln;
}


template<typename Scalar, typename StorageIndex>
void  myIncompleteLUT2<Scalar,StorageIndex>::free()
{
  if(nRow == 0)
    return;
  
  int phase = -1, error = 0;
  double ddum;

  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &nRow, array, &csr[0], &col[0], &perm[0], &NRHS,
           IPARM, &MSGLVL, &ddum, &ddum, &error, DPARM);
  
  array = NULL;

  currentStatus = EMPTY;

  return;
}
    

}


#endif








