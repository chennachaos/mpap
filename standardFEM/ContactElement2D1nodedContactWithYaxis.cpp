

#include "ContactElement2D1nodedContactWithYaxis.h"
#include "BasisFunctionsLagrange.h"
#include "GeomDataLagrange.h"
#include "SolutionData.h"

using namespace std;


ContactElement2D1nodedContactWithYaxis::ContactElement2D1nodedContactWithYaxis()
{
  ndof   = 1;
  ndim   = 2;
  npElem = 1;
  nlbf   = npElem;
  nsize  = npElem*ndof;
}


ContactElement2D1nodedContactWithYaxis::~ContactElement2D1nodedContactWithYaxis()
{
}

void ContactElement2D1nodedContactWithYaxis::prepareElemData()
{
  return;
}


int ContactElement2D1nodedContactWithYaxis::calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal)
{
    if(Klocal.rows() != nsize)
    {
      Klocal.resize(nsize, nsize);
      Flocal.resize(nsize);
    }
    Klocal.setZero();
    Flocal.setZero();


    int ii, jj, ind;
    double Yorig, disp, lamn;
    double  y1, y2, fact, tol=1.e-12, gn, af, cn=1.0, g0, d1;

    Yorig = GeomData->NodePosOrig[nodeNums[0]][1];
    disp  =  SolnData->var1Cur[nodeNums[0]];
    lamn  =  SolnData->var1Cur[nodeNums[0]];

    int type=1;

    if(type == 1)  // scenario #1
    {
      y1 = 0.00575;
      g0 = 0.0001;
    }
    else // relief valve
    {
      y1 = 21.1; // top-most point on the lower block
      g0 = 0.05; // initial gap
    }

    y2 = y1 + g0 + disp;

    gn = - disp; // gn is penetration variable

    //cout << " displacement = " << disp << '\t' << lamn << endl;
    //cout << " penetration = " << gn << '\t' << lamn << endl;

    af = SolnData->td(2);
    d1 = SolnData->td(5);
  fact = 1.0/SolnData->td(10);

  if( (lamn + cn*gn) > tol)
  {
    Klocal(0,1) -= af*cn;
    Klocal(1,0) -= af*fact;

    Flocal(0)   += lamn;
    Flocal(1)   -= gn*cn;  // residual: contact force

  }
  else
  {
    Klocal(1,1) -= af;
    Flocal(1)   += lamn;
  }

  //printMatrix(Klocal);

  return 1;
}




/*
void ContactElement2D1nodedContactWithYaxis::assembleMatrixAndVector(int start, SparseMatrixXd& mtx, double* rhs)
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


