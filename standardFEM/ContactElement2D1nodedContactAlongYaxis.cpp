

#include "ContactElement2D1nodedContactAlongYaxis.h"
#include "BasisFunctionsLagrange.h"
#include "GeomDataLagrange.h"
#include "SolutionData.h"

using namespace std;


ContactElement2D1nodedContactAlongYaxis::ContactElement2D1nodedContactAlongYaxis()
{
  ndof   = 1;
  ndim   = 2;
  npElem = 1;
  nlbf   = npElem;
  nsize  = 2; // to account for the Lagrange multiplier
}


ContactElement2D1nodedContactAlongYaxis::~ContactElement2D1nodedContactAlongYaxis()
{
}

void ContactElement2D1nodedContactAlongYaxis::prepareElemData()
{
  return;
}

int ContactElement2D1nodedContactAlongYaxis::calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal, bool firstIter)
{
    if(Klocal.rows() != nsize)
    {
      Klocal.resize(nsize, nsize);
      Flocal.resize(nsize);
    }
    Klocal.setZero();
    Flocal.setZero();

    double  tol   = 1.e-12, cn=1.0;
    double  Yorig =  GeomData->NodePosOrig[nodeNums[0]][1];
    double  disp  =  SolnData->var1Cur[nodeNums[0]*2+1];
    double  lamn  =  SolnData->var1Cur[GeomData->assy4r[forAssyVec[1]]];

    //cout << forAssyVec[0] << '\t' << forAssyVec[1] << endl;
    //cout << GeomData->assy4r[forAssyVec[1]] << endl;
    //cout << endl;

    double  gn = -disp; // gn is penetration variable

    double  af = SolnData->td(2);
    double  d1 = SolnData->td(5);

    if( (lamn + cn*gn) > tol)
    {
      //cout << " displacement = " << disp << '\t' << lamn << endl;

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

    return 0;
}


/*
int ContactElement2D1nodedContactAlongYaxis::calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal, bool firstIter)
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
*/


