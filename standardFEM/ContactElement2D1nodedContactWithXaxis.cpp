
#include "ContactElement2D1nodedContactWithXaxis.h"
#include "BasisFunctionsLagrange.h"

#include "ImmersedSolid.h"


using namespace std;


ContactElement2D1nodedContactWithXaxis::ContactElement2D1nodedContactWithXaxis()
{
  ndof   = 1;
  ndim   = 2;
  npElem = 1;
  nlbf   = npElem;
  nsize  = 2; // to account for the Lagrange multiplier
}


ContactElement2D1nodedContactWithXaxis::~ContactElement2D1nodedContactWithXaxis()
{
}


void ContactElement2D1nodedContactWithXaxis::prepareElemData()
{
  return;
}


int ContactElement2D1nodedContactWithXaxis::calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal)
{
    if(Klocal.rows() != nsize)
    {
      Klocal.resize(nsize, nsize);
      Flocal.resize(nsize);
    }
    Klocal.setZero();
    Flocal.setZero();


    int ii, jj, ind;
    double  Yorig, disp, lamn;
    double  y1, y2, fact, tol=1.e-12, gn, af, cn=1.0, g0, d1;

    Yorig =  GeomData->NodePosOrig[nodeNums[0]][1];
    disp  =  SolnData->var1Cur[nodeNums[0]*2-1];
    lamn  =  SolnData->var1Cur[GeomData->assy4r[forAssyVec[1]]];

    gn = - disp; // gn is penetration variable

    //cout << " displacement = " << disp << '\t' << lamn << endl;

    af = SolnData->td(2);
    d1 = SolnData->td(5);

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

    //printMatrix(Klocal);

  return 1;
}




/*
void ContactElement2D1nodedContactWithXaxis::AssembleMatrixAndVector(int start, SparseMatrixXd& mtx, double* rhs)
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


