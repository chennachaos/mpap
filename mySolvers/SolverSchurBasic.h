#ifndef incl_SolverSchurBasic_h
#define incl_SolverSchurBasic_h

#include "headersBasic.h"

using namespace std;
using namespace Eigen;

typedef  SparseMatrix<double, RowMajor>  SparseMatrixXd;

int SolverSchurCG(const SparseMatrixXd& Amat, const SparseMatrixXd& Bmat, const SparseMatrixXd& Cmat, const VectorXd& f3, 
VectorXd& var1, VectorXd& var2, int niter, double TOL, bool usePreCond)
{

}


#endif
