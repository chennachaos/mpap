
#include "LagrangeElem1Dbar2Node.h"

#include "Debug.h"
#include "MpapTime.h"
#include "ComputerTime.h"
#include "GeomDataLagrange.h"
#include "SolutionData.h"

using namespace std;

extern ComputerTime       computerTime;
extern MpapTime mpapTime;


LagrangeElem1Dbar2Node::LagrangeElem1Dbar2Node()
{
  degree = 1;
  ndof   = 1;
  ndim   = 2;
  npElem = 2;
  nlbf = npElem;
  nsize  = npElem*ndof;

  if (debug) cout << " constructor LagrangeElem1Dbar2Node\n\n";
}



LagrangeElem1Dbar2Node::~LagrangeElem1Dbar2Node()
{
  if (debug) cout << " destructor LagrangeElem1Dbar2Node\n\n";
}



void LagrangeElem1Dbar2Node::prepareElemData()
{
  LagrangeElement::prepareElemData();
  
  //Klocal.resize(nsize, nsize);
  //Flocal.resize(nsize);

  return;
}




void LagrangeElem1Dbar2Node::prepareElemData2()
{
  int ii, jj, ind, kk, ind2;
  
  if(!(SolnData->STAGGERED) )
  {
    //cout << " aaaaaaaaaaaaaa " << npElem*ndim << endl;

    ind = npElem*ndim;

    forAssyVec2.resize(ind);

    kk = 0;
    for(ii=0;ii<npElem;ii++)
    {
      ind2 = nodeNums[ii]*ndim;
      for(jj=0;jj<ndim;jj++)
      {
        //cout << ii << '\t' << jj << '\t' << kk << '\t' << ind2 << endl;
        forAssyVec2[kk++] = ind2 + jj;
      }
    }
  }
  
  return;
}



int LagrangeElem1Dbar2Node::calcLoadVector()
{
  return 0;
}



int LagrangeElem1Dbar2Node::calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal, bool firstIter)
{
    int  ii, jj, gp, kk, ll, TI, TIp1, TIp2, count, TJ, TJp1, TJp2;
    
    double  bf, af, am, dt, d1, rho, aa;
    double  A, E, EA, dvol, fact, fact1, fact2, fact3, fact4;
    double  x2, y2, x1, y1, u1, u2, w1, w2, ddu[2], ddw[2], ddbt[2], h;

    VectorXd  dispC(2), accC(2), velC(2);
    MatrixXd  Mlocal(2,2);

    //elmDat = &(SolidSolnData->ElemProp[elmType].data[0]);
    //matDat = &(SolidSolnData->MatlProp[matType].data[0]);

    //bf   = elmDat[4];
    //E    = elmDat[5];
    //A    = elmDat[6];
    bf = 1.0;
    E  = 1.0;
    A  = 1.0;

    //cout << " material constants... "  << EA << '\t' << EI << '\t' << GA << endl;
    
    dt = SolnData->td(0);
    af = SolnData->td(2);
    af = 1.0;
    am = SolnData->td(1);
    d1 = SolnData->td(5);
    aa = SolnData->td(10);
    
    //  rotate nodal displacements and compute nodal positions on element axis
    
    //cout << " ccccccccccc " << endl;
    //cout << nodeNums[0] << '\t' << nodeNums[1] << endl;
    //cout << SolidSolnData->NodePosOrig[nodeNums[0]][0] << '\t' << SolidSolnData->NodePosOrig[nodeNums[0]][1] << endl;

    x1 = GeomData->NodePosOrig[nodeNums[0]][0];
    y1 = GeomData->NodePosOrig[nodeNums[0]][1];
    x2 = GeomData->NodePosOrig[nodeNums[1]][0];
    y2 = GeomData->NodePosOrig[nodeNums[1]][1];
    
    //cout << " ccccccccccc " << endl;

    ii = nodeNums[0]*ndof;
    jj = nodeNums[1]*ndof;

    //dispC(0) =  SolidSolnData->dispCur[nodeNums[0]];
    //dispC(1) =  SolidSolnData->dispCur[nodeNums[1]];

    //velC(0)  =  SolidSolnData->veloCur[nodeNums[0]];
    //velC(1)  =  SolidSolnData->veloCur[nodeNums[1]];

    //accC(0)  =  SolidSolnData->acceCur[nodeNums[0]];
    //accC(1)  =  SolidSolnData->acceCur[nodeNums[1]];

    dispC(0) =  SolnData->var1Cur[nodeNums[0]];
    dispC(1) =  SolnData->var1Cur[nodeNums[1]];

    // compute the orientation of the element

    h = x2-x1;

    //cout << x1 << '\t' << x2 << '\t' << h << endl;

    Klocal.setZero();
    Flocal.setZero();
    
    fact = E*A/h;

    Klocal(0,0) =  fact;
    Klocal(0,1) = -fact;
    Klocal(1,0) = -fact;
    Klocal(1,1) =  fact;

    //printMatrix(Klocal);	printf("\n\n\n");
    //printVector(dispC); printf("\n\n");

    //body forces
    
    Flocal(0) = 0.5*h*bf;
    Flocal(1) = 0.5*h*bf;

    Flocal -= Klocal*dispC;
    Klocal = af*Klocal;

    //printMatrix(Klocal); printf("\n\n");
    //printVector(Flocal); printf("\n\n");

   return 0;
}
//





int LagrangeElem1Dbar2Node::calcInternalForces()
{
  return 0;
}



void LagrangeElem1Dbar2Node::discreteContourplot(int vartype, int varindex, int index, int nCol, double umin, double umax)
{
  return;
}


void LagrangeElem1Dbar2Node::projectToKnots(bool extrapolateFlag, int vartype, int varindex, int index)
{
  return;
}


void LagrangeElem1Dbar2Node::projectStress(int varindex, double* outval)
{
  return;
}


void LagrangeElem1Dbar2Node::projectStrain(int vartype, int varindex, double* outval)
{
  return;
}



void LagrangeElem1Dbar2Node::projectIntVar(int index, double* outval)
{
   return;
}


int LagrangeElem1Dbar2Node::calcOutput(double u1, double v1)
{
  return 0;
}


void LagrangeElem1Dbar2Node::toPostprocess(int vartype, int varindex, int type, SparseMatrixXd&  coeffMat, VectorXd& rhsVec)
{
  return;
}

