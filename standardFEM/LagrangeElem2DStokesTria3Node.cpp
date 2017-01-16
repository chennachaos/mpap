
#include "LagrangeElem2DStokesTria3Node.h"

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


LagrangeElem2DStokesTria3Node::LagrangeElem2DStokesTria3Node()
{
  degree = 1;
  npElem = 3;
  nlbf   = 3;
  ndof   = 3;
  nsize  = nlbf*ndof;

  if (debug) cout << " constructor LagrangeElem2DStokesTria3Node\n\n";
}

LagrangeElem2DStokesTria3Node::~LagrangeElem2DStokesTria3Node()
{
  if (debug) cout << " destructor LagrangeElem2DStokesTria3Node\n\n";
}


void LagrangeElem2DStokesTria3Node::prepareElemData()
{
  LagrangeElement::prepareElemData();

//  Klocal.resize(nsize, nsize);
//  Flocal.resize(nsize);

  return;
}


void LagrangeElem2DStokesTria3Node::prepareElemData2()
{
  return;
}



void LagrangeElem2DStokesTria3Node::resetMatrixAndVector()
{
  return;
}



int LagrangeElem2DStokesTria3Node::calcLoadVector()
{
  return 0;
}



int LagrangeElem2DStokesTria3Node::calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal)
{
//   char fct[] = "LagrangeElem2DNavierStokes4Node::calcStiffnessAndResidual";
//   computerTime.go(fct);

    //Stokes2DEx1 analy;
    //Kovasznay  analy;
    //analy.SetPressure(0.0);

    int ii, jj, gp1, gp2, TI, TIp1, TIp2, count, TJ, TJp1, TJp2;

    double  uu, vv, Jac, dvol, b1, b2, b3, b4, b5, b6, b7, b8, xx, yy, acceFact, HH, dist, rho, mu;
    double  pres, Da, Db, af, am, d1, c1, muTaf, rad, urdr, urdr2, tau[3], volume;
    double  fact, fact1, fact2, bb1, bb2, param[2], h2, stabParam, *elmDat, x[3], y[3];

    for(ii=0;ii<npElem;ii++)
    {
      x[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      y[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    VectorXd  N(nlbf), dN_dx(nlbf), dN_dy(nlbf);

    VectorXd  res(3), dp(2), Du(2), vel(2), velDot(2), force(2), res2(2), rStab(2);
    MatrixXd  grad(2,2), stress(2,2), Dj(2,3);

    elmDat = &(SolnData->ElemProp[elmType].data[0]);
    //matDat = &(SolnData->MatlProp[matType].data[0]);

    rho = elmDat[4];
    mu  = elmDat[5];

    bool axsy = false;
    axsy = ((int)elmDat[2] == 1);

    af = SolnData->td(2);
    am = SolnData->td(3);
    acceFact = rho*am*SolnData->td(9);
    muTaf = mu*af;

    volume = 0.5*(x[0]*(y[1]-y[2]) + x[1]*(y[2]-y[0]) + x[2]*(y[0]-y[1]));

    //cout << " volume = " << volume << endl;
    h2 = 4.0*volume/PI;

    stabParam = h2/(4.0*mu);
    tau[0] = elmDat[8]*stabParam;  // SUPG
    tau[1] = elmDat[9]*stabParam;  // PSPG
    tau[2] = elmDat[10]*stabParam; // LSIC

    int nGP1 = 2;
    int nGP2 = 2;

    getGaussPointsTriangle(nGP1, GeomData->gausspoints1, GeomData->gausspoints2, GeomData->gaussweights1);

    if(Klocal.rows() != nsize)
    {
      Klocal.resize(nsize, nsize);
      Flocal.resize(nsize);
    }
    Klocal.setZero();
    Flocal.setZero();
    
    //cout << " AAAAAAAAAA " << endl;
    //cout << nGP1 << '\t' << nGP2 << endl;
    //cout << " AAAAAAAAAA " << endl;
    //cout << nGP1 << '\t' << nGP2 << endl;

    for(gp1=0;gp1<nGP1;gp1++)
    {
          param[0] = GeomData->gausspoints1[gp1];
          param[1] = GeomData->gausspoints2[gp1];

          GeomData->computeBasisFunctions2D(0, 1, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

          dvol = GeomData->gaussweights1[gp1] * Jac;
          
          //printf("%12.6f \t %12.6f \t %12.6f \n", volume, h, Jac);
          //printf("%12.6f \t %12.6f \t %12.6f \t %12.6f \n", matB(0,0), matB(0,1), matB(1,0), matB(1,1));
          //printf("%12.6f \t %12.6f \t %12.6f \t %12.6f \n", matG(0,0), matG(0,1), matG(1,0), matG(1,1));

          vel(0) = computeValueCur(0, N);
          vel(1) = computeValueCur(1, N);

          grad(0,0) = computeValueCur(0, dN_dx);
          grad(0,1) = computeValueCur(0, dN_dy);
          grad(1,0) = computeValueCur(1, dN_dx);
          grad(1,1) = computeValueCur(1, dN_dy);

          Du.setZero();

          pres   = computeValue(2, N);
          dp(0)  = computeValue(2, dN_dx);
          dp(1)  = computeValue(2, dN_dy);

          velDot(0) = computeValueDotCur(0, N);
          velDot(1) = computeValueDotCur(1, N);
          //velDot.setZero();

          // this is pseudo-stress
          stress = mu*grad;
          stress(0,0) -= pres;
          stress(1,1) -= pres;

          //cout << xx << '\t' << yy << endl;
          force.setZero();
          //force(0) = analy.computeForce(0, xx, yy);
          //force(1) = analy.computeForce(1, xx, yy);
          //cout << force(0) << '\t' << force(1) << endl;

          //force(0) = computeForce(0, N);
          //force(1) = computeForce(1, N);

          res2(0) = rho*(velDot(0) - force(0)) ;
          res2(1) = rho*(velDot(1) - force(1)) ;

          rStab(0) = res2(0) - mu*Du(0) + dp(0) ;
          rStab(1) = res2(1) - mu*Du(1) + dp(1) ;

          if(axsy)
          {
            xx = computeGeomOrig(0, N);
            yy = computeGeomOrig(1, N);

            rad = xx;

            urdr  = vel(0)/rad;
            urdr2 = urdr/rad;
            dvol *= (2.0*PI*rad);

            rStab(0) -= mu*(grad(0,0)/rad - urdr2 );
            rStab(1) -= mu*(grad(1,0)/rad );
          }

          tau[0] *= elmDat[8];  // SUPG
          tau[1] *= elmDat[9];  // PSPG
          tau[2] *= elmDat[10]; // LSIC

          for(ii=0;ii<nlbf;ii++)
          {
            TI   = ndof*ii;
            TIp1 = TI+1;
            TIp2 = TI+2;

            b1 = dN_dx[ii]*dvol;
            b2 = dN_dy[ii]*dvol;
            b4 = N[ii]*dvol;

            b5 = muTaf*b1;
            b6 = muTaf*b2;
            b8 = af*b4;

            for(jj=0;jj<nlbf;jj++)
            {
              TJ   = ndof*jj;
              TJp1 = TJ+1;
              TJp2 = TJ+2;

              fact2 = rho*acceFact*N(jj);

              // time acceleration term
              fact = b4*fact2 ;

              // diffusion term
              fact += ( b5*dN_dx(jj)+b6*dN_dy(jj) );

              Klocal(TI,   TJ)   += fact;
              Klocal(TIp1, TJp1) += fact;

              // pressure term
              Klocal(TI,   TJp2) -= (b1*N(jj));
              Klocal(TIp1, TJp2) -= (b2*N(jj));

              // continuity equation
              Klocal(TIp2, TJ)   += (b8*dN_dx(jj));
              Klocal(TIp2, TJp1) += (b8*dN_dy(jj));

              // PSPG stabilisation terms

              Dj(0,0) = fact2;
              Dj(0,1) = 0.0;
              Dj(0,2) = dN_dx(jj);
              Dj(1,0) = 0.0;
              Dj(1,1) = fact2;
              Dj(1,2) = dN_dy(jj);

              if(axsy)
              {
                Dj(0,0) -= muTaf*(dN_dx(jj)/rad - N(jj)/rad/rad);
                Dj(1,1) -= muTaf*(dN_dx(jj)/rad);
              }

              // PSPG stabilisation
              Klocal(TIp2, TJ)   += (b1*Dj(0,0) + b2*Dj(1,0))*tau[1];
              Klocal(TIp2, TJp1) += (b1*Dj(0,1) + b2*Dj(1,1))*tau[1];
              Klocal(TIp2, TJp2) += (b1*Dj(0,2) + b2*Dj(1,2))*tau[1];

              // LSIC stabilisation

              fact = af*rho*tau[2];

              Klocal(TI,   TJ)   += (b1*fact*dN_dx(jj));
              Klocal(TI,   TJp1) += (b1*fact*dN_dy(jj));

              Klocal(TIp1, TJ)   += (b2*fact*dN_dx(jj));
              Klocal(TIp1, TJp1) += (b2*fact*dN_dy(jj));

              if(axsy)
              {
                  // diffusion term
                  Klocal(TI, TJ)     += (b4 * (mu/rad/rad) * (af*N(jj)) );
                  Klocal(TI, TJp2)   -= (b4 * N(jj)/rad);

                  // continuity equation
                  Klocal(TIp2, TJ)   += (b4 * af*N(jj)/rad);
              }
            }

            Flocal(TI)   -= (b4*res2(0) + b1*stress(0,0) + b2*stress(0,1) );
            Flocal(TIp1) -= (b4*res2(1) + b1*stress(1,0) + b2*stress(1,1) );
            Flocal(TIp2) -= (b4*grad.trace());

            // PSPG stabilisation terms
            Flocal(TIp2) -= (tau[1]*(b1*rStab(0)+b2*rStab(1)));

            // LSIC stabilisation terms
            fact2 = tau[2]*rho*grad.trace();

            Flocal(TI)   -= b1*fact2;
            Flocal(TIp1) -= b2*fact2;

            if(axsy)
            {
                Flocal(TI)   -= (b4 * (mu/rad/rad) * vel(0) );
                Flocal(TI)   += (b4 * pres/rad);
                Flocal(TIp2) -= (b4 * vel(0)/rad);
            }
          }
  }//gp1
  
  //printMatrix(Klocal); printf("\n\n\n");
  //printVector(Flocal); printf("\n\n\n");


  //cout << nGP1 << '\t' << nGP2 << endl;

  return 0;
}






int LagrangeElem2DStokesTria3Node::calcInternalForces()
{
  return 0;
}



void LagrangeElem2DStokesTria3Node::discreteContourplot(int vartype, int varindex, int index, int nCol, double umin, double umax)
{
  return;
}


void LagrangeElem2DStokesTria3Node::projectToKnots(bool extrapolateFlag, int vartype, int varindex, int index)
{
  return;
}


void LagrangeElem2DStokesTria3Node::projectStress(int varindex, double* outval)
{
  return;
}



void LagrangeElem2DStokesTria3Node::projectStrain(int vartype, int varindex, double* outval)
{
  return;
}



void LagrangeElem2DStokesTria3Node::projectIntVar(int index, double* outval)
{
  return;
}


int LagrangeElem2DStokesTria3Node::calcOutput(double u1, double v1)
{
  return 0;
}



void LagrangeElem2DStokesTria3Node::toPostprocess(int vartype, int varindex, int type, SparseMatrixXd&  coeffMat, VectorXd& rhsVec)
{
  return;
}


