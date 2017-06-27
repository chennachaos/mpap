
#include "LagrangeElem3DStokesTet4Node.h"

#include "Debug.h"
#include "MpapTime.h"
#include "ComputerTime.h"
#include "GeomDataLagrange.h"
#include "SolutionData.h"
#include "Functions.h"
#include "QuadratureUtil.h"
#include "myConstants.h"

using namespace std;

extern ComputerTime       computerTime;
extern MpapTime mpapTime;


LagrangeElem3DStokesTet4Node::LagrangeElem3DStokesTet4Node()
{
  ndim   = 3;
  degree = 1;
  npElem = 4;
  nlbf   = 4;
  ndof   = 4;
  nsize  = nlbf*ndof;

  if (debug) cout << " constructor LagrangeElem3DStokesTet4Node\n\n";
}

LagrangeElem3DStokesTet4Node::~LagrangeElem3DStokesTet4Node()
{
  if (debug) cout << " destructor LagrangeElem3DStokesTet4Node\n\n";
}


void LagrangeElem3DStokesTet4Node::prepareElemData()
{
  LagrangeElement::prepareElemData();

  return;
}


void LagrangeElem3DStokesTet4Node::prepareElemData2()
{
  return;
}


int LagrangeElem3DStokesTet4Node::calcLoadVector()
{
  return 0;
}



int LagrangeElem3DStokesTet4Node::calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal)
{
  double  dvol, b1, b2, b3, b4, b5, b6, b7, b8;
  double  Jac, totvol, pres, fact, fact2, *elmDat;
  double  param[3], bforce[3], tau[3];

  VectorXd  N(nlbf), dN_dx(nlbf), dN_dy(nlbf), dN_dz(nlbf);
  VectorXd  vel(3), velDot(3), res2(3), Du(3), dp(3), rStab(3), force(3);
  MatrixXd  grad(3,3), gradN(3,3), stress(3,3),  Dj(3,4), forVolume(4,4);

  int   index, ii, jj, gp, TI, TIp1, TIp2, TIp3, TJ, TJp1, TJp2, TJp3;

  elmDat = &(SolnData->ElemProp[elmType].data[0]);

  double rho = elmDat[3] ;
  double mu  = elmDat[4] ;

  bforce[0]  = 0.0;
  bforce[1]  = 0.0;
  bforce[2]  = 0.0;

  //bforce[0]  = elmDat[6]*timeFunction[0].prop ;
  //bforce[1]  = elmDat[7]*timeFunction[0].prop ;
  //bforce[2]  = elmDat[7]*timeFunction[0].prop ;

  double dt = mpapTime.dt;
  double af = SolnData->td(2);
  double am = SolnData->td(3);
  double acceFact = am*SolnData->td(9);
  double muTaf = mu*af;

  //printf(" %14.12f \t %14.12f \t %14.12f \t %14.12f \n ", BULK, mu, bforce[0], bforce[1]);
  //printf(" %14.12f \t  %14.12f \t %14.12f \t %14.12f \t %14.12f \n ", rho, mu, af, am, acceFact);
  //printf(" %5d \t %5d \t %5d \t %5d \n ", nodeNums[0], nodeNums[1], nodeNums[2], nodeNums[3]);

  double xNode[4], yNode[4], zNode[4], xx, yy, zz;

  for(ii=0;ii<npElem;ii++)
  {
    xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
    yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    zNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][2];
  }

  //forVolume(0,0) = 1.0; forVolume(0,1) = xNode[0]; forVolume(0,2) = yNode[0]; forVolume(0,3) = zNode[0];
  //forVolume(1,0) = 1.0; forVolume(1,1) = xNode[1]; forVolume(1,2) = yNode[1]; forVolume(1,3) = zNode[1];
  //forVolume(2,0) = 1.0; forVolume(2,1) = xNode[2]; forVolume(2,2) = yNode[2]; forVolume(2,3) = zNode[2];
  //forVolume(3,0) = 1.0; forVolume(3,1) = xNode[3]; forVolume(3,2) = yNode[3]; forVolume(3,3) = zNode[3];

  forVolume(0,0) = 1.0;      forVolume(0,1) = 1.0;      forVolume(0,2) = 1.0;      forVolume(0,3) = 1.0;
  forVolume(1,0) = xNode[0]; forVolume(1,1) = xNode[1]; forVolume(1,2) = xNode[2]; forVolume(1,3) = xNode[3];
  forVolume(2,0) = yNode[0]; forVolume(2,1) = yNode[1]; forVolume(2,2) = yNode[2]; forVolume(2,3) = yNode[3];
  forVolume(3,0) = zNode[0]; forVolume(3,1) = zNode[1]; forVolume(3,2) = zNode[2]; forVolume(3,3) = zNode[3];

  totvol = forVolume.determinant()/6.0;

  //if(elenum < 100)
    //cout << " volume = " << totvol << endl;

  double h = pow(6.0*totvol/PI, 1.0/3.0);
  
  double stabParam = h*h/(4.0*mu);
  tau[0] = elmDat[8]*stabParam;  // SUPG
  tau[1] = elmDat[9]*stabParam;  // PSPG
  tau[2] = elmDat[10]*stabParam; // LSIC

  //cout << " totvol = " << totvol << '\t' << tau[1] << endl;

  //printf(" %14.12f \t %14.12f \t %14.12f \t %14.12f \n ", xNode[0], xNode[1], xNode[2], xNode[3]);
  //printf(" %14.12f \t %14.12f \t %14.12f \t %14.12f \n\n\n ", yNode[0], yNode[1], yNode[2], yNode[3]);

  if(Klocal.rows() != nsize)
  {
    Klocal.resize(nsize, nsize);
    Flocal.resize(nsize);
  }
  Klocal.setZero();
  Flocal.setZero();

  //cout << " finite = " << finite << endl;

  //nGP = (int) elmDat[0] ;
  nGP = 1;

  vector<double>  gausspoints1, gausspoints2, gausspoints3, gaussweights;

  getGaussPointsTet(nGP, gausspoints1, gausspoints2, gausspoints3, gaussweights);

  for(gp=0;gp<nGP;gp++)
  {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];
        param[2] = gausspoints3[gp];

        GeomData->computeBasisFunctions3D(0, 1, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);

        dvol = Jac * gaussweights[gp];

          vel(0) = computeValueCur(0, N);
          vel(1) = computeValueCur(1, N);
          vel(2) = computeValueCur(2, N);

          grad(0,0) = computeValueCur(0, dN_dx);
          grad(0,1) = computeValueCur(0, dN_dy);
          grad(0,2) = computeValueCur(0, dN_dz);

          grad(1,0) = computeValueCur(1, dN_dx);
          grad(1,1) = computeValueCur(1, dN_dy);
          grad(1,2) = computeValueCur(1, dN_dz);

          grad(2,0) = computeValueCur(2, dN_dx);
          grad(2,1) = computeValueCur(2, dN_dy);
          grad(2,2) = computeValueCur(2, dN_dz);

          Du.setZero();

          pres   = computeValue(3, N);
          dp(0)  = computeValue(3, dN_dx);
          dp(1)  = computeValue(3, dN_dy);
          dp(2)  = computeValue(3, dN_dz);

          velDot(0) = computeValueDotCur(0, N);
          velDot(1) = computeValueDotCur(1, N);
          velDot(2) = computeValueDotCur(2, N);


          // this is pseudo-stress
          stress = mu*grad;
          stress(0,0) -= pres;
          stress(1,1) -= pres;
          stress(2,2) -= pres;

          //cout << xx << '\t' << yy << endl;
          force.setZero();
          //force(0) = analy.computeForce(0, xx, yy);
          //force(1) = analy.computeForce(1, xx, yy);
          //cout << force(0) << '\t' << force(1) << endl;

          res2(0) = rho*(velDot(0) - force(0)) ;
          res2(1) = rho*(velDot(1) - force(1)) ;
          res2(2) = rho*(velDot(2) - force(2)) ;

          rStab(0) = res2(0) - mu*Du(0) + dp(0) ;
          rStab(1) = res2(1) - mu*Du(1) + dp(1) ;
          rStab(2) = res2(2) - mu*Du(2) + dp(2) ;


          for(ii=0;ii<nlbf;ii++)
          {
            TI   = 4*ii;
            TIp1 = TI+1;
            TIp2 = TI+2;
            TIp3 = TI+3;

            b1 = dN_dx[ii]*dvol;
            b2 = dN_dy[ii]*dvol;
            b3 = dN_dz[ii]*dvol;
            b4 = N[ii]*dvol;

            b5 = muTaf*b1;
            b6 = muTaf*b2;
            b7 = muTaf*b3;
            b8 = af*b4;

            for(jj=0;jj<nlbf;jj++)
            {
              TJ   = 4*jj;
              TJp1 = TJ+1;
              TJp2 = TJ+2;
              TJp3 = TJ+3;

              fact2 = rho*acceFact*N(jj);
              //cout << " fact2 = " << fact2 << endl;

              // time acceleration term
              fact = b4*fact2 ;

              // diffusion term
              fact += ( b5*dN_dx(jj)+b6*dN_dy(jj)+b7*dN_dz(jj) );

              Klocal(TI,   TJ)   += fact;
              Klocal(TIp1, TJp1) += fact;
              Klocal(TIp2, TJp2) += fact;

              // pressure term
              Klocal(TI,   TJp3) -= (b1*N(jj));
              Klocal(TIp1, TJp3) -= (b2*N(jj));
              Klocal(TIp2, TJp3) -= (b3*N(jj));

              // continuity equation
              Klocal(TIp3, TJ)   += (b8*dN_dx(jj));
              Klocal(TIp3, TJp1) += (b8*dN_dy(jj));
              Klocal(TIp3, TJp2) += (b8*dN_dz(jj));

              // PSPG stabilisation terms

              Dj(0,0) = fact2;
              Dj(0,1) = 0.0;
              Dj(0,2) = 0.0;
              Dj(0,3) = af*dN_dx(jj);

              Dj(1,0) = 0.0;
              Dj(1,1) = fact2;
              Dj(1,2) = 0.0;
              Dj(1,3) = af*dN_dy(jj);

              Dj(2,0) = 0.0;
              Dj(2,1) = 0.0;
              Dj(2,2) = fact2;
              Dj(2,3) = af*dN_dz(jj);


              // PSPG stabilisation
              Klocal(TIp3, TJ)   += (b1*Dj(0,0) + b2*Dj(1,0) + b3*Dj(2,0))*tau[1];
              Klocal(TIp3, TJp1) += (b1*Dj(0,1) + b2*Dj(1,1) + b3*Dj(2,1))*tau[1];
              Klocal(TIp3, TJp2) += (b1*Dj(0,2) + b2*Dj(1,2) + b3*Dj(2,2))*tau[1];
              Klocal(TIp3, TJp3) += (b1*Dj(0,3) + b2*Dj(1,3) + b3*Dj(2,3))*tau[1];

              // LSIC stabilisation

              fact = af*rho*tau[2];

              Klocal(TI,   TJ)   += (b1*fact*dN_dx(jj));
              Klocal(TI,   TJp1) += (b1*fact*dN_dy(jj));
              Klocal(TI,   TJp2) += (b1*fact*dN_dz(jj));

              Klocal(TIp1, TJ)   += (b2*fact*dN_dx(jj));
              Klocal(TIp1, TJp1) += (b2*fact*dN_dy(jj));
              Klocal(TIp1, TJp2) += (b2*fact*dN_dz(jj));

              Klocal(TIp2, TJ)   += (b3*fact*dN_dx(jj));
              Klocal(TIp2, TJp1) += (b3*fact*dN_dy(jj));
              Klocal(TIp2, TJp2) += (b3*fact*dN_dz(jj));
            }

            Flocal(TI)   -= (b4*res2(0) + b1*stress(0,0) + b2*stress(0,1) + b3*stress(0,2) );
            Flocal(TIp1) -= (b4*res2(1) + b1*stress(1,0) + b2*stress(1,1) + b3*stress(1,2) );
            Flocal(TIp2) -= (b4*res2(2) + b1*stress(2,0) + b2*stress(2,1) + b3*stress(2,2) );
            Flocal(TIp3) -= (b4*grad.trace());

            // PSPG stabilisation terms
            Flocal(TIp3) -= (tau[1]*(b1*rStab(0)+b2*rStab(1)+b3*rStab(2)));

            // LSIC stabilisation terms
            fact2 = tau[2]*rho*grad.trace();

            Flocal(TI)   -= b1*fact2;
            Flocal(TIp1) -= b2*fact2;
            Flocal(TIp2) -= b3*fact2;
          }
  }//gp

  //printMatrix(Klocal); printf("\n\n\n");
  //printVector(Flocal); printf("\n\n\n");

  return 0;
}






int LagrangeElem3DStokesTet4Node::calcInternalForces()
{
  return 0;
}



void LagrangeElem3DStokesTet4Node::discreteContourplot(int vartype, int varindex, int index, int nCol, double umin, double umax)
{
  return;
}


void LagrangeElem3DStokesTet4Node::projectToKnots(bool extrapolateFlag, int vartype, int varindex, int index)
{
  return;
}


void LagrangeElem3DStokesTet4Node::projectStress(int varindex, double* outval)
{
  return;
}



void LagrangeElem3DStokesTet4Node::projectStrain(int vartype, int varindex, double* outval)
{
  return;
}



void LagrangeElem3DStokesTet4Node::projectIntVar(int index, double* outval)
{
  return;
}


int LagrangeElem3DStokesTet4Node::calcOutput(double u1, double v1)
{
  return 0;
}



void LagrangeElem3DStokesTet4Node::toPostprocess(int vartype, int varindex, int type, SparseMatrixXd&  coeffMat, VectorXd& rhsVec)
{
  return;
}


