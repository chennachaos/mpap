
#include "LagrangeElem2DPoissonTria3Node.h"

#include "Debug.h"
#include "MpapTime.h"
#include "ComputerTime.h"
#include "GeomDataLagrange.h"
#include "SolutionData.h"
#include "QuadratureUtil.h"

using namespace std;

extern ComputerTime       computerTime;
extern MpapTime mpapTime;


LagrangeElem2DPoissonTria3Node::LagrangeElem2DPoissonTria3Node()
{
  degree = 1;
  npElem = 3;
  nlbf   = 3;
  ndof   = 1;
  nsize  = 3;

  if (debug) cout << " constructor LagrangeElem2DPoissonTria3Node\n\n";
}

LagrangeElem2DPoissonTria3Node::~LagrangeElem2DPoissonTria3Node()
{
  if (debug) cout << " destructor LagrangeElem2DPoissonTria3Node\n\n";
}


void LagrangeElem2DPoissonTria3Node::prepareElemData()
{
  LagrangeElement::prepareElemData();

  return;
}


void LagrangeElem2DPoissonTria3Node::prepareElemData2()
{
  return;
}



int LagrangeElem2DPoissonTria3Node::calcLoadVector()
{
  return 0;
}




/*
int LagrangeElem2DPoissonTria3Node::calcStiffnessAndResidual(double* xco, double* yco, double* SolnData)
{
  double  area, x1, x2, x3, y1, y2, y3, y23, y31, y12, x32, x13, x21, kx, ky;
  
  kx = 1.0;
  ky = 1.0;

  MatrixXd  Bmat(3,2);
  VectorXd  dispC(3);

  dispC.setZero();

    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

  x1 = GeomData->NodePosOrig[nodeNums[0]][0];
  y1 = GeomData->NodePosOrig[nodeNums[0]][1];

  x2 = GeomData->NodePosOrig[nodeNums[1]][0];
  y2 = GeomData->NodePosOrig[nodeNums[1]][1];

  x3 = GeomData->NodePosOrig[nodeNums[2]][0];
  y3 = GeomData->NodePosOrig[nodeNums[2]][1];

  x1 = xco[0];
  x2 = xco[1];
  x3 = xco[2];
  y1 = yco[0];
  y2 = yco[1];
  y3 = yco[2];

//    cout << x1 << '\t' << y1 << endl;
//    cout << x2 << '\t' << y2 << endl;
//    cout << x3 << '\t' << y3 << endl;

  solnElem(0) =  SolnData->var1Cur[nodeNums[0]];
  solnElem(1) =  SolnData->var1Cur[nodeNums[1]];
  solnElem(2) =  SolnData->var1Cur[nodeNums[2]];

  dispC(0) =  SolnData[nodeNums[0]];
  dispC(1) =  SolnData[nodeNums[1]];
  dispC(2) =  SolnData[nodeNums[2]];
  
  area = 0.5*(x2*y3 - x3*y2 + x3*y1 - x1*y3 + x1*y2 - x2*y1);
  
  //cout << " area = " << area << endl;

  x13 = x1 - x3;  x21 = x2 - x1;  x32 = x3 - x2;
  y31 = y3 - y1;  y12 = y1 - y2;  y23 = y2 - y3;

  Klocal.setZero();
  Flocal.setZero();
  
  Bmat(0,0) = y23; Bmat(1,0) = y31; Bmat(2,0) = y12;
  Bmat(0,1) = x32; Bmat(1,1) = x13; Bmat(2,1) = x21;
  
  Bmat /= (2.0*area);
  
  Klocal = (Bmat*(kx*area)) * Bmat.transpose();
  Flocal -= Klocal*dispC;
  
//if(subdomId == 1)
//{
  //printMatrix(Klocal);   printf("\n\n\n");
  //printVector(Flocal);   printf("\n\n\n");
//}

//  Klocal(0,0) = 0.5;  Klocal(0,1) = 0.0;  Klocal(0,2) = -0.5;
//  Klocal(1,0) = 0.0;  Klocal(1,1) = 0.5;  Klocal(1,2) = -0.5;
//  Klocal(2,0) =-0.5;  Klocal(2,1) =-0.5;  Klocal(2,2) =  1.0;

  forAssyVec = nodeNums;

  return 0;
}
*/


//int LagrangeElem2DPoissonTria3Node::calcStiffnessAndResidual(double* xco, double* yco, double* SolnData)

int LagrangeElem2DPoissonTria3Node::calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal)
{
  int ii, jj, gp;

  double  area, x1, x2, x3, y1, y2, y3, y23, y31, y12, x32, x13, x21;
  double  du[2], param[2], b1, b2, b3, Jac, dvol, fact, force, xx, yy, xNode[3], yNode[3];
  
  VectorXd  N(nlbf), dN_dx(nlbf), dN_dy(nlbf);
  VectorXd  solnElem(3);
  solnElem.setZero();

  for(ii=0;ii<npElem;ii++)
  {
    xNode[ii] = GeomData->NodePosOrig[SolnData->node_map_new_to_old[nodeNums[ii]]][0];
    yNode[ii] = GeomData->NodePosOrig[SolnData->node_map_new_to_old[nodeNums[ii]]][1];
  }

  double totVolume = 0.5*(xNode[0]*(yNode[1]-yNode[2]) + xNode[1]*(yNode[2]-yNode[0]) + xNode[2]*(yNode[0]-yNode[1]));

  //cout << " volume = " << totVolume << endl;

  solnElem(0) =  SolnData->var1Cur[nodeNums[0]];
  solnElem(1) =  SolnData->var1Cur[nodeNums[1]];
  solnElem(2) =  SolnData->var1Cur[nodeNums[2]];

  nGP = 1;

    vector<double>  gausspoints1, gausspoints2, gaussweights;

    getGaussPointsTriangle(nGP, gausspoints1, gausspoints2, gaussweights);

    if(Klocal.rows() != nsize)
    {
      Klocal.resize(nsize, nsize);
      Flocal.resize(nsize);
    }

    Klocal.setZero();
    Flocal.setZero();

    for(gp=0;gp<nGP;gp++)
    {
          param[0] = gausspoints1[gp];
          param[1] = gausspoints2[gp];

          //computeBasisFunctions2D(0, 1, degree, param, nodeNums, xNode, yNode, &N(0), &dN_dx(0), &dN_dy(0), Jac);
          GeomData->computeBasisFunctions2D(0, 1, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

          dvol = gaussweights[gp] * Jac;
          
          //printf("%12.6f \t %12.6f \t %12.6f \n", Jac, dvol, totVolume);
          //printf("%12.6f \t %12.6f \t %12.6f \t %12.6f \n", matB(0,0), matB(0,1), matB(1,0), matB(1,1));
          //printf("%12.6f \t %12.6f \t %12.6f \t %12.6f \n", matG(0,0), matG(0,1), matG(1,0), matG(1,1));

          xx = yy = 0.0;
          du[0] = du[1] = 0.0;

          for(ii=0;ii<nlbf;ii++)
          {
            xx += xNode[ii]*N[ii];
            yy += yNode[ii]*N[ii];

            du[0] += solnElem(ii)*dN_dx(ii);
            du[1] += solnElem(ii)*dN_dy(ii);
          }

          //cout << xx << '\t' << yy << endl;
          //force = -cos(PI*xx)*cos(PI*yy);
          force = 0.0;

          for(ii=0;ii<nlbf;ii++)
          {
            b1 = dN_dx[ii]*dvol;
            b2 = dN_dy[ii]*dvol;
            b3 = N[ii]*dvol;

            Flocal(ii)  += (b3*force - b1*du[0] - b2*du[1] );

            for(jj=0;jj<nlbf;jj++)
            {
              Klocal(ii, jj)   += (b1*dN_dx[jj] + b2*dN_dy[jj]);
            }
          }
  }//gp1
  
  //printMatrix(Klocal); printf("\n\n\n");
  //printVector(Flocal); printf("\n\n\n");

  return 1;
}




int LagrangeElem2DPoissonTria3Node::calcInternalForces()
{
  return 0;
}



void LagrangeElem2DPoissonTria3Node::discreteContourplot(int vartype, int varindex, int index, int nCol, double umin, double umax)
{
  return;
}


void LagrangeElem2DPoissonTria3Node::projectToKnots(bool extrapolateFlag, int vartype, int varindex, int index)
{
  return;
}


void LagrangeElem2DPoissonTria3Node::projectStress(int varindex, double* outval)
{
  return;
}



void LagrangeElem2DPoissonTria3Node::projectStrain(int vartype, int varindex, double* outval)
{
  return;
}



void LagrangeElem2DPoissonTria3Node::projectIntVar(int index, double* outval)
{
  return;
}


int LagrangeElem2DPoissonTria3Node::calcOutput(double u1, double v1)
{
  return 0;
}



void LagrangeElem2DPoissonTria3Node::toPostprocess(int vartype, int varindex, int type, SparseMatrixXd&  coeffMat, VectorXd& rhsVec)
{
  return;
}


