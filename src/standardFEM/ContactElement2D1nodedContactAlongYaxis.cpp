

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
    //cout << GeomData->assy4r[forAssyVec[1]] << endl;     cout << endl;

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
      //cout << " af = " << af << '\t' << lamn << endl;

      Klocal(1,1) += af;
      Flocal(1)   -= lamn;
    }

    //printMatrix(Klocal);

    return 0;
}

