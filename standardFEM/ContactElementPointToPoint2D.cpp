
#include "ContactElementPointToPoint2D.h"
#include "BasisFunctionsLagrange.h"

#include "ImmersedSolid.h"


using namespace std;


ContactElementPointToPoint2D::ContactElementPointToPoint2D()
{
  DIM = 2;
}


ContactElementPointToPoint2D::~ContactElementPointToPoint2D()
{
}

void ContactElementPointToPoint2D::prepareElemData()
{
  int ii, jj;

  DIM = 2;
  nNode = 2;
  ndof = 3;

  size = nNode*ndof;

  nodePos.resize(nNode);
  boundaryConds.resize(nNode);

  for(ii=0;ii<nNode;ii++)
  {
    nodePos[ii].resize(DIM);
    boundaryConds[ii].resize(ndof);
  }

  nodePos[0][0] = 0.0;
  nodePos[0][1] = 0.00575;

  nodePos[1][0] = 0.0;
  nodePos[1][1] = 0.005752;

  boundaryConds[0][0] = -1; // x-disp
  boundaryConds[0][1] = -1; // y-disp
  boundaryConds[0][2] = -1; // x-lambda

  boundaryConds[1][0] = -1; // x-lambda
  boundaryConds[1][1] =  1; // x-lambda
  boundaryConds[1][2] = -1; // x-lambda

  
  totalDOF=0;
  for(ii=0;ii<nNode;ii++)
  {
    for(jj=0;jj<DIM;jj++)
    {
      if(boundaryConds[ii][jj] == 1)
        assy4r.push_back(totalDOF++);
    }
  }
  
  size = 2;

  Klocal.resize(size,size);
  Flocal.resize(size);

  return;
}


void ContactElementPointToPoint2D::initialiseDOFvalues()
{
  return;
}


void ContactElementPointToPoint2D::reset()
{
  return;
}



void ContactElementPointToPoint2D::calcStiffnessAndResidual(int ind1, int ind2, double inp1, double inp2)
{
  int ii, jj, ind;
  double  y1, y2, fact, tol=1.e-12, disp, lamn, gn, af, cn=1.0, g0;

  Klocal.setZero();
  Flocal.setZero();

  disp = SolidSolnData->var1Cur[1];
  //lamn = SolnData.var1Cur[1];
  lamn = SolidSolnData->var1DotCur[3];
  //lamn = FluidSolnData->var4Cur(1);
  
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

  af   = SolidSolnData->td[2];
  fact = 1.0/SolidSolnData->td[10];

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

  return;
}






void ContactElementPointToPoint2D::assembleElementVector(int ind, bool flag, double* rhs)
{
  return;
}


void ContactElementPointToPoint2D::assembleElementMatrix(int index, SparseMatrixXd& mtx)
{
  return;
}


void ContactElementPointToPoint2D::assembleMatrixAndVector(int start, SparseMatrixXd& mtx, double* rhs)
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



