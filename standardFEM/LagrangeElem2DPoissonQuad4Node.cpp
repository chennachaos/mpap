
#include "LagrangeElem2DPoissonQuad4Node.h"

#include "Debug.h"
#include "MpapTime.h"
#include "ComputerTime.h"
#include "GeomDataLagrange.h"
#include "SolutionData.h"

using namespace std;

extern ComputerTime       computerTime;
extern MpapTime mpapTime;


LagrangeElem2DPoissonQuad4Node::LagrangeElem2DPoissonQuad4Node()
{
  degree = 1;
  npElem = 4;
  nlbf   = 4;
  ndof   = 1;
  nsize  = 4;

  if (debug) cout << " constructor LagrangeElem2DPoissonQuad4Node\n\n";
}

LagrangeElem2DPoissonQuad4Node::~LagrangeElem2DPoissonQuad4Node()
{
  if (debug) cout << " destructor LagrangeElem2DPoissonQuad4Node\n\n";
}


void LagrangeElem2DPoissonQuad4Node::prepareElemData()
{
  LagrangeElement::prepareElemData();

  return;
}


void LagrangeElem2DPoissonQuad4Node::prepareElemData2()
{
  return;
}


int LagrangeElem2DPoissonQuad4Node::calcLoadVector()
{
  return 0;
}



int LagrangeElem2DPoissonQuad4Node::calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal)
{
//   char fct[] = "LagrangeElem2DPoissonQuad4Node::calcStiffnessAndResidual";
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

  int nGP1 = 2;
  int nGP2 = 2;

  for(gp2=0;gp2<nGP2;gp2++)
  {
    param[1] = GeomData->gausspoints2[gp2];
  for(gp1=0;gp1<nGP1;gp1++)
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






int LagrangeElem2DPoissonQuad4Node::calcInternalForces()
{
  return 0;
}



void LagrangeElem2DPoissonQuad4Node::discreteContourplot(int vartype, int varindex, int index, int nCol, double umin, double umax)
{
  return;
}


void LagrangeElem2DPoissonQuad4Node::projectToKnots(bool extrapolateFlag, int vartype, int varindex, int index)
{
  return;
}


void LagrangeElem2DPoissonQuad4Node::projectStress(int varindex, double* outval)
{
  return;
}



void LagrangeElem2DPoissonQuad4Node::projectStrain(int vartype, int varindex, double* outval)
{
  return;
}



void LagrangeElem2DPoissonQuad4Node::projectIntVar(int index, double* outval)
{
  return;
}


int LagrangeElem2DPoissonQuad4Node::calcOutput(double u1, double v1)
{
  return 0;
}



void LagrangeElem2DPoissonQuad4Node::toPostprocess(int vartype, int varindex, int type, SparseMatrixXd&  coeffMat, VectorXd& rhsVec)
{
  return;
}


