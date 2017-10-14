

#ifndef EIGEN_MY_INCOMPLETE_LUT_H
#define EIGEN_MY_INCOMPLETE_LUT_H

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

  typedef int Index;


/** \internal
  * Compute a quick-sort split of a vector 
  * On output, the vector row is permuted such that its elements satisfy
  * abs(row(i)) >= abs(row(ncut)) if i<ncut
  * abs(row(i)) <= abs(row(ncut)) if i>ncut 
  * \param row The vector of values
  * \param ind The array of index for the elements in @p row
  * \param ncut  The number of largest elements to keep
  **/ 
template <typename VectorV, typename VectorI>
Index QuickSplit2(VectorV &row, VectorI &ind, Index ncut)
{
  typedef typename VectorV::RealScalar RealScalar;
  using std::swap;
  using std::abs;
  Index mid;
  Index n = row.size(); /* length of the vector */
  Index first, last, j ;

  RealScalar abskey;

  ncut--; /* to fit the zero-based indices */
  first = 0; 
  last = n-1; 
  if(ncut < first || ncut > last )
    return 0;

  do
  {
    mid = first; 
    abskey = abs(row(mid)); 
    for(j = first + 1; j <= last; j++)
    {
      if( abs(row(j)) > abskey)
      {
        ++mid;
        swap(row(mid), row(j));
        swap(ind(mid), ind(j));
      }
    }
    /* Interchange for the pivot element */
    swap(row(mid), row(first));
    swap(ind(mid), ind(first));
    
    if(mid > ncut)
      last = mid - 1;
    else if (mid < ncut )
      first = mid + 1; 

  }while (mid != ncut );
  
  return 0; /* mid is equal to ncut */ 
}




template <typename _Scalar, typename _StorageIndex = int>
class  myIncompleteLUT : public SparseSolverBase<myIncompleteLUT<_Scalar, _StorageIndex> >
{
  protected:
    typedef SparseSolverBase<myIncompleteLUT> Base;
    using Base::m_isInitialized;
  public:
    typedef _Scalar Scalar;
    typedef _StorageIndex StorageIndex;
    typedef typename NumTraits<Scalar>::Real RealScalar;
    typedef Matrix<Scalar,Dynamic,1> Vector;
    typedef Matrix<StorageIndex,Dynamic,1> VectorI;
    typedef SparseMatrix<Scalar,RowMajor,StorageIndex> FactorType;
    //typedef SparseMatrix<Scalar,RowMajor> FactorType;

    StorageIndex   nRow, nCol;
    int  currentStatus;
    VectorXd  soln;

  public:
    
    // this typedef is only to export the scalar type and compile-time dimensions to solve_retval
    typedef Matrix<Scalar,Dynamic,Dynamic> MatrixType;
    
    myIncompleteLUT()
      : m_droptol(NumTraits<Scalar>::dummy_precision()), m_fillfactor(10),
        m_analysisIsOk(false), m_factorizationIsOk(false), nRow(0)
    {}

    //~myIncompleteLUT();

    void free();

    template<typename MatrixType>
    explicit  myIncompleteLUT(const MatrixType& mat, const RealScalar& droptol=NumTraits<Scalar>::dummy_precision(), int fillfactor = 10)
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
    int factorize(const MatrixType& amat);
    
    /**
      * Compute an incomplete LU factorization with dual threshold on the matrix mat
      * No pivoting is done in this version
      * 
      **/
    template<typename MatrixType>
    myIncompleteLUT&  compute(const MatrixType& amat)
    {
      analyzePattern(amat);
      factorize(amat);
      
      return *this;
    }

    void setDroptol(const RealScalar& droptol); 
    void setFillfactor(int fillfactor); 

    /*
    VectorXd&  solve(const VectorXd&  rhsVec)
    {
      soln = m_Pinv * rhsVec;
      soln = m_lu.template triangularView<UnitLower>().solve(soln);
      soln = m_lu.template triangularView<Upper>().solve(soln);
      soln = m_P * soln; 
    }
    */

    template<typename Rhs, typename Dest>
    void _solve_impl(const Rhs& b, Dest& x) const
    {
      x = m_Pinv * b;
      x = m_lu.template triangularView<UnitLower>().solve(x);
      x = m_lu.template triangularView<Upper>().solve(x);
      x = m_P * x; 

      //x = m_lu.template triangularView<UnitLower>().solve(b);
      //x = m_lu.template triangularView<Upper>().solve(x);
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
 * Set control parameter droptol
 *  \param droptol   Drop any element whose magnitude is less than this tolerance 
 **/ 
template<typename Scalar, typename StorageIndex>
void myIncompleteLUT<Scalar,StorageIndex>::setDroptol(const RealScalar& droptol)
{
  this->m_droptol = droptol;   
}



/**
 * Set control parameter fillfactor
 * \param fillfactor  This is used to compute the  number @p fill_in of largest elements to keep on each row. 
 **/ 
template<typename Scalar, typename StorageIndex>
void myIncompleteLUT<Scalar,StorageIndex>::setFillfactor(int fillfactor)
{
  this->m_fillfactor = fillfactor;   
}



template <typename Scalar, typename StorageIndex>
template<typename _MatrixType>
void myIncompleteLUT<Scalar,StorageIndex>::analyzePattern(const _MatrixType& amat)
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

  // Then compute the AMD ordering...
  
  //internal::minimum_degree_ordering<Scalar, StorageIndex>(amat.selfadjointView<Lower>(), m_P);
  //internal::minimum_degree_ordering<Scalar, StorageIndex>(amat, m_P);

  //AMDOrdering<StorageIndex> ordering;
  
  //ordering(mat1, m_P);
  
  //m_P.setIdentity();
  
  //std::cout << m_P.indices() << endl;

  m_Pinv  = m_P.inverse(); // ... and the inverse permutation

  m_analysisIsOk = true;
  m_factorizationIsOk = false;
  m_isInitialized = false;
  
  return;
}



template <typename Scalar, typename StorageIndex>
template<typename _MatrixType>
int  myIncompleteLUT<Scalar,StorageIndex>::factorize(const _MatrixType& amat)
{
//
  using std::sqrt;
  using std::max;
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
  
  cout << amat.nonZeros() << '\t' << mat.nonZeros() << endl;
  
  const int* rowPtr   = amat.outerIndexPtr();
  const int* colIndx  = amat.innerIndexPtr();
  const double* array = amat.valuePtr();

  const int* rowPtr2   = amat.outerIndexPtr();
  const int* colIndx2  = amat.innerIndexPtr();
  const double* array2 = amat.valuePtr();

  /*
  FILE * pFile;
  char name [100];

  pFile = fopen ("Poisson2.dat","w");

  fprintf(pFile, "%9d \t %9d \t %9d \n", n, n, 0 );

  for(int k=0; k<mat.outerSize(); ++k)
  for(SparseMatrixXd::InnerIterator it(mat,k); it; ++it)
    fprintf(pFile, "%9d \t %9d \t %20.16f \n", it.row(), it.col(), it.value());

  fclose (pFile);
  */

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
  Index ii;
  m_lu.reserve(n * (nnzL + nnzU + 1));

/*
  int max1 = -1000000;
  int max2 = -1000000;
  
  Index ikk;

  for(ii = 0; ii < n; ii++)
  {
    //ikk =  rowPtr[ii+1]-rowPtr[ii];
    //for(Index ik=0; ik<ikk; ik++)
      max1 = max(max1, (rowPtr[ii+1]-rowPtr[ii]) );
      max2 = max(max2, (rowPtr2[ii+1]-rowPtr2[ii]) );
  }
  
  cout << max1 << '\t' << max2 << endl;
*/


  // global loop over the rows of the sparse matrix
  //#pragma omp parallel for default(shared) private(ii)
  for(ii = 0; ii < n; ii++)
  {
    //Vector u(n) ;     // real values of the row -- maximum size is n --
    //VectorI ju(n);   // column position of the values in u -- maximum size  is n
    //VectorI jr(n);   // Indicate the position of the nonzero elements in the vector u 
    //                 //  -- A zero location is indicated by -1

    // Initialization
    //jr.fill(-1);
    //ju.fill(0);
    //u.fill(0);


    // 1 - copy the lower and the upper part of the row i of mat in the working vector u

    Index sizeu = 1; // number of nonzero elements in the upper part of the current row
    Index sizel = 0; // number of nonzero elements in the lower part of the current row

    ju(ii)    = ii;
    u(ii)     = 0;
    jr(ii)    = ii;

    RealScalar rownorm = 0.0;

    // Iterate through the current row ii

    Index  k1 = rowPtr[ii];
    Index  ikk = (rowPtr[ii+1]-rowPtr[ii]);
    Index  ik;
    
    //cout << k1 << '\t' << ikk << endl;
    
    //#pragma omp threadprivate(sizel, sizeu)  
    //#pragma omp parallel  default(shared) reduction(+: rownorm, sizel, sizeu) private(ik) //shared(ii, k1, ikk, ju, u, jr, rowPtr, colIndx, array)
    //{
    //#pragma omp for 
    //for(Index ik=OuterStarts[ii]; ik<OuterStarts[ii+1]; ik++)
    for(ik=0; ik<ikk; ik++)
    {
      Index k2 = k1+ik;
      Index k = colIndx[k2];
      RealScalar  temp = array[k2];

      //#pragma omp critical
        //cout << ii << '\t' << ik << '\t' << k << '\t' << sizel << '\t' << sizeu << endl;

      if (k < ii)
      {
        // copy the lower part
        ju(sizel) = k;
        u(sizel) = temp;
        jr(k) = sizel;
        //#pragma omp atomic
          ++sizel;
        //sizel += 1;
      }
      else if (k == ii)
      {
        u(ii) = temp;
      }
      else
      {
        // copy the upper part
        Index jpos = ii + sizeu;
        ju(jpos) = k;
        u(jpos) = temp;
        jr(k) = jpos;
        //#pragma omp atomic
          ++sizeu;
        //sizeu += 1;
      }

      rownorm += numext::abs2(temp);

    }
    //}

    // 2 - detect possible zero row
    if(rownorm==0)
    {
      m_info = NumericalIssue;
      cout << " myIncompleteLUT .... Numerical issues ... rownorm=0 " << endl;

      return -1;
    }
    // Take the 2-norm of the current row as a relative tolerance
    rownorm = sqrt(rownorm);

    // 3 - eliminate the previous nonzero rows
    Index jj = 0;
    Index len = 0;
    while(jj < sizel)
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
        jr(minrow) = jj;
        jr(j) = k;
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
          ju(newpos) = j;
          u(newpos) = -prod;
          jr(j) = newpos;
        }
        else
          u(jpos) -= prod;
      }
      // store the pivot element
      u(len)  = fact;
      ju(len) = minrow;
      ++len;

      jj++;
    } // end of the elimination on the row ii

    // reset the upper part of the pointer jr to zero
    //#pragma omp parallel for
    for(Index k=0; k<sizeu; k++)
      jr(ju(ii+k)) = -1;

    // 4 - partially sort and insert the elements in the m_lu matrix

    // sort the L-part of the row
    sizel = len;
    len = (std::min)(sizel, nnzL);

    typename Vector::SegmentReturnType ul(u.segment(0, sizel));

    typename VectorI::SegmentReturnType jul(ju.segment(0, sizel));

    internal::QuickSplit(ul, jul, len);

    // store the largest m_fill elements of the L part
    //#pragma omp critical
    //{
    m_lu.startVec(ii);
    for(Index k = 0; k < len; k++)
      m_lu.insertBackByOuterInnerUnordered(ii,ju(k)) = u(k);
    //}

    // store the diagonal element
    // apply a shifting rule to avoid zero pivots (we are doing an incomplete factorization)
    if (u(ii) == Scalar(0))
      u(ii) = sqrt(m_droptol) * rownorm;

    m_lu.insertBackByOuterInnerUnordered(ii, ii) = u(ii);

    // sort the U-part of the row
    // apply the dropping rule first
    len = 0;
    for(Index k=1; k<sizeu; k++)
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
    //#pragma omp critical
    for(Index k= ii+1; k < ii+len; k++)
      m_lu.insertBackByOuterInnerUnordered(ii,ju(k)) = u(k);
  }

  m_lu.finalize();
  m_lu.makeCompressed();
  
  cout <<  m_lu.nonZeros()  << endl;

  m_factorizationIsOk = true;
  m_isInitialized = true;
  m_info = Success;

  return 1;
}



template<typename Scalar, typename StorageIndex>
void  myIncompleteLUT<Scalar,StorageIndex>::free()
{
  return;
}


}



#endif








