
#include "LagrangeElem3DPoissonHex8Node.h"

#include "Debug.h"
#include "MpapTime.h"
#include "ComputerTime.h"
#include "GeomDataLagrange.h"
#include "SolutionData.h"

using namespace std;

extern ComputerTime       computerTime;
extern MpapTime mpapTime;


LagrangeElem3DPoissonHex8Node::LagrangeElem3DPoissonHex8Node()
{
  degree = 1;
  npElem = 4;
  nlbf   = 4;
  ndof   = 1;
  nsize  = 4;

  if (debug) cout << " constructor LagrangeElem3DPoissonHex8Node\n\n";
}

LagrangeElem3DPoissonHex8Node::~LagrangeElem3DPoissonHex8Node()
{
  if (debug) cout << " destructor LagrangeElem3DPoissonHex8Node\n\n";
}


void LagrangeElem3DPoissonHex8Node::prepareElemData()
{
  LagrangeElement::prepareElemData();

  return;
}


void LagrangeElem3DPoissonHex8Node::prepareElemData2()
{
  return;
}


int LagrangeElem3DPoissonHex8Node::calcLoadVector()
{
  return 0;
}



int LagrangeElem3DPoissonHex8Node::calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal)
{
//   char fct[] = "LagrangeElem3DPoissonHex8Node::calcStiffnessAndResidual";
//   computerTime.go(fct);

  int gp1, gp2;
  double  fact, dvol, Jac, param[2];

  VectorXd  N(nlbf), dN_dx(nlbf), dN_dy(nlbf);

    if(Klocal.rows() != nsize)
    {
      Klocal.resize(nsize, nsize);
      Flocal.resize(nsize);
    }
    Klocal.setZero();
    Flocal.setZero();


  nGP = 2;

  for(gp2=0; gp2<nGP; gp2++)
  {
    param[1] = GeomData->gausspoints2[gp2];
  for(gp1=0; gp1<nGP; gp1++)
  {
        param[0] = GeomData->gausspoints1[gp1];

        GeomData->computeBasisFunctions2D(0, 2, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

        dvol = GeomData->gaussweights2[gp2] * GeomData->gaussweights1[gp1] * Jac;

        //cout << " AAAAAAAAAAA " << endl;
        Flocal += (dvol*(N*fact - dN_dx*computeValueCur(0, dN_dx) - dN_dy*computeValueCur(0, dN_dy)));
        //cout << " AAAAAAAAAAA " << endl;
        Klocal += (dvol*(dN_dx*dN_dx.transpose()+dN_dy*dN_dy.transpose()));
  }//gp1
  }//gp2
  
  //printMatrix(Klocal);   printf("\n\n\n");
  //printVector(Flocal);   printf("\n\n\n");

  return 0;
}






int LagrangeElem3DPoissonHex8Node::calcInternalForces()
{
  return 0;
}



void LagrangeElem3DPoissonHex8Node::discreteContourplot(int vartype, int varindex, int index, int nCol, double umin, double umax)
{
  return;
}


void LagrangeElem3DPoissonHex8Node::projectToKnots(bool extrapolateFlag, int vartype, int varindex, int index)
{
  return;
}


void LagrangeElem3DPoissonHex8Node::projectStress(int varindex, double* outval)
{
  return;
}



void LagrangeElem3DPoissonHex8Node::projectStrain(int vartype, int varindex, double* outval)
{
  return;
}



void LagrangeElem3DPoissonHex8Node::projectIntVar(int index, double* outval)
{
  return;
}


int LagrangeElem3DPoissonHex8Node::calcOutput(double u1, double v1)
{
  return 0;
}



void LagrangeElem3DPoissonHex8Node::toPostprocess(int vartype, int varindex, int type, SparseMatrixXd&  coeffMat, VectorXd& rhsVec)
{
  return;
}


