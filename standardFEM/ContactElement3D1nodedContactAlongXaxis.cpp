

#include "ContactElement3D1nodedContactAlongXaxis.h"
#include "BasisFunctionsLagrange.h"
#include "GeomDataLagrange.h"
#include "SolutionData.h"

using namespace std;


ContactElement3D1nodedContactAlongXaxis::ContactElement3D1nodedContactAlongXaxis()
{
  ndof   = 1;
  ndim   = 3;
  npElem = 1;
  nlbf   = npElem;
  nsize  = 2; // to account for the Lagrange multiplier
}


ContactElement3D1nodedContactAlongXaxis::~ContactElement3D1nodedContactAlongXaxis()
{
}

void ContactElement3D1nodedContactAlongXaxis::prepareElemData()
{
  return;
}


int ContactElement3D1nodedContactAlongXaxis::calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal, bool firstIter)
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
    double  disp  =  SolnData->var1Cur[nodeNums[0]*3+0];
    double  lamn  =  SolnData->var1Cur[GeomData->assy4r[forAssyVec[1]]];

    double  gn = 0.0 - disp; // gn is penetration variable

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

    //printMatrix(Klocal);

    return 1;
}


