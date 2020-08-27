
#include "Debug.h"
#include "SolverEigen.h"
#include "SolverTime.h"
#include "ComputerTime.h"
#include "util.h"


//#include <Eigen/SuperLUSupport>
#include <Eigen/SparseExtra>
#include <Eigen/IterativeSolvers>


extern SolverTime      solverTime;
extern ComputerTime    computerTime;

using namespace std;
using namespace Eigen;



void SolverEigen::computeConditionNumber()
{
/*
    MatrixXd  globalK;
    
    globalK.resize(nRow, nCol);
    globalK.setZero();

    int k, ii, jj;

    for(k=0; k<mtx.outerSize(); ++k)
    {
      for(SparseMatrixXd::InnerIterator it(mtx,k); it; ++it)
      {
        ii = it.row();
        jj = it.col();
        
        //cout << ii << '\t' << jj << '\t' << it.value() << endl;

        globalK.coeffRef(ii, jj) = it.value();
      }
    }


    VectorXd sing_vals = globalK.jacobiSvd().singularValues();

    //printf("\n Matrix condition number = %12.6E \n", sing_vals(0)/sing_vals(sing_vals.size()-1) );
    
    //printf("\n Minimum eigenvalue = %12.6f \n", sing_vals.minCoeff() );
    //printf("\n Minimum eigenvalue = %12.6f \n", sing_vals.maxCoeff() );
    printf("\n Matrix condition number = %12.6E \n", sing_vals.maxCoeff() / sing_vals.minCoeff() );
    printf("\n\n\n\n");
*/

  //myCondNumMatlab(mtx);

  return;
}







