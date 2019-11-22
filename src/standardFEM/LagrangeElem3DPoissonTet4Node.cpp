
#include "LagrangeElem3DPoissonTet4Node.h"

#include "Debug.h"
#include "MpapTime.h"
#include "ComputerTime.h"
#include "GeomDataLagrange.h"
#include "SolutionData.h"
#include "QuadratureUtil.h"
#include "Functions.h"


using namespace std;

extern ComputerTime       computerTime;
extern MpapTime mpapTime;


LagrangeElem3DPoissonTet4Node::LagrangeElem3DPoissonTet4Node()
{
  degree = 1;
  npElem = 4;
  nlbf   = 4;
  ndof   = 1;
  nsize  = 4;

  if (debug) cout << " constructor LagrangeElem3DPoissonTet4Node\n\n";
}

LagrangeElem3DPoissonTet4Node::~LagrangeElem3DPoissonTet4Node()
{
  if (debug) cout << " destructor LagrangeElem3DPoissonTet4Node\n\n";
}


void LagrangeElem3DPoissonTet4Node::prepareElemData()
{
  LagrangeElement::prepareElemData();

  return;
}


void LagrangeElem3DPoissonTet4Node::prepareElemData2()
{
  return;
}


int LagrangeElem3DPoissonTet4Node::calcLoadVector()
{
  return 0;
}



int LagrangeElem3DPoissonTet4Node::calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal, bool firstIter)
{
//   char fct[] = "LagrangeElem3DPoissonTet4Node::calcStiffnessAndResidual";
//   computerTime.go(fct);

  PoissonEx2  analy;

  int gp, ii;
  double  force, dvol, Jac, param[3], xx, yy, zz, xNode[4], yNode[4], zNode[4];

  for(ii=0;ii<npElem;ii++)
  {
    xNode[ii] = GeomData->NodePosOrig[SolnData->node_map_new_to_old[nodeNums[ii]]][0];
    yNode[ii] = GeomData->NodePosOrig[SolnData->node_map_new_to_old[nodeNums[ii]]][1];
    zNode[ii] = GeomData->NodePosOrig[SolnData->node_map_new_to_old[nodeNums[ii]]][2];
  }


  VectorXd  N(nlbf), dN_dx(nlbf), dN_dy(nlbf), dN_dz(nlbf);

    if(Klocal.rows() != nsize)
    {
      Klocal.resize(nsize, nsize);
      Flocal.resize(nsize);
    }
    Klocal.setZero();
    Flocal.setZero();

  elmDat = &(SolnData->ElemProp[elmType].data[0]);

  nGP = (int) elmDat[0] ;

  vector<double>  gausspoints1, gausspoints2, gausspoints3, gaussweights;

  getGaussPointsTet(nGP, gausspoints1, gausspoints2, gausspoints3, gaussweights);


  for(gp=0;gp<nGP;gp++)
  {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];
        param[2] = gausspoints3[gp];

        GeomData->computeBasisFunctions3D(0, 1, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);

        dvol = gaussweights[gp] * Jac;
        
          xx = yy = zz = 0.0;
          //du[0] = du[1] = 0.0;

          for(ii=0;ii<nlbf;ii++)
          {
            xx += xNode[ii]*N[ii];
            yy += yNode[ii]*N[ii];
            zz += zNode[ii]*N[ii];

            //du[0] += solnElem(ii)*dN_dx(ii);
            //du[1] += solnElem(ii)*dN_dy(ii);
          }

        //force = analy.computeForce(0, xx, yy, zz);
        force = 0.0;
        //cout << " AAAAAAAAAAA " << endl;
        Flocal += (dvol*(N*force - dN_dx*computeValueCur(0, dN_dx) - dN_dy*computeValueCur(0, dN_dy) - dN_dz*computeValueCur(0, dN_dz)));
        //cout << " AAAAAAAAAAA " << endl;
        Klocal += (dvol*(dN_dx*dN_dx.transpose()+dN_dy*dN_dy.transpose()+dN_dz*dN_dz.transpose()));
  }//gp
  
  //printMatrix(Klocal);   printf("\n\n\n");
  //printVector(Flocal);   printf("\n\n\n");

  return 0;
}






int LagrangeElem3DPoissonTet4Node::calcInternalForces()
{
  return 0;
}



void LagrangeElem3DPoissonTet4Node::discreteContourplot(int vartype, int varindex, int index, int nCol, double umin, double umax)
{
  return;
}


void LagrangeElem3DPoissonTet4Node::projectToKnots(bool extrapolateFlag, int vartype, int varindex, int index)
{
  return;
}


void LagrangeElem3DPoissonTet4Node::projectStress(int varindex, double* outval)
{
  return;
}



void LagrangeElem3DPoissonTet4Node::projectStrain(int vartype, int varindex, double* outval)
{
  return;
}



void LagrangeElem3DPoissonTet4Node::projectIntVar(int index, double* outval)
{
  return;
}


int LagrangeElem3DPoissonTet4Node::calcOutput(double u1, double v1)
{
  return 0;
}



void LagrangeElem3DPoissonTet4Node::toPostprocess(int vartype, int varindex, int type, SparseMatrixXd&  coeffMat, VectorXd& rhsVec)
{
  return;
}


