
#include "ContactElement2D1nodedContactAlongXaxis.h"
#include "BasisFunctionsLagrange.h"
#include "GeomDataLagrange.h"
#include "SolutionData.h"


using namespace std;


ContactElement2D1nodedContactAlongXaxis::ContactElement2D1nodedContactAlongXaxis()
{
  ndof   = 1;
  ndim   = 2;
  npElem = 1;
  nlbf   = npElem;
  nsize  = 2; // to account for the Lagrange multiplier
}


ContactElement2D1nodedContactAlongXaxis::~ContactElement2D1nodedContactAlongXaxis()
{
}


void ContactElement2D1nodedContactAlongXaxis::prepareElemData()
{
  return;
}


int ContactElement2D1nodedContactAlongXaxis::calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal, bool firstIter)
{
    if(Klocal.rows() != nsize)
    {
      Klocal.resize(nsize, nsize);
      Flocal.resize(nsize);
    }
    Klocal.setZero();
    Flocal.setZero();


    double  tol=1.e-12, cn=1.0;

    double  Xorig =  GeomData->NodePosOrig[nodeNums[0]][0];
    double  disp  =  SolnData->var1Cur[nodeNums[0]*2+0];
    double  lamn  =  SolnData->var1Cur[GeomData->assy4r[forAssyVec[1]]];

    double  gn = 0.0-disp; // gn is penetration variable

    //cout << " displacement = " << disp << '\t' << lamn << endl;

    double  af = SolnData->td(2);
    double  d1 = SolnData->td(5);


    if( (lamn + cn*gn) > tol)
    {
      Klocal(0,1) -= af;
      Klocal(1,0) -= af;

      Flocal(0)   += lamn;
      Flocal(1)   -= gn;  // residual: contact force
    }
    else
    {
      Klocal(1,1) += af;
      Flocal(1)   -= lamn;
    }

    //printMatrix(Klocal);  printf("\n\n\n");  printVector(Flocal);

    return 1;
}




/*
void ContactElement2D1nodedContactAlongXaxis::assembleMatrixAndVector(int start, SparseMatrixXd& mtx, double* rhs)
{
  int ii, jj, r;

  //cout << " start " << start << endl;

  for(ii=0;ii<size;ii++)
  {
    r = start+ii;

    rhs[r] += Flocal[ii];

    for(jj=0;jj<size;jj++)
    {
      mtx.coeffRef(r, start+jj) += Klocal(ii,jj);
    }
  }
  
  return;
}
*/


