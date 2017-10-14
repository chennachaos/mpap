
#include "LagrangeElem2DStokesQuad4Node.h"

#include "Debug.h"
#include "MpapTime.h"
#include "ComputerTime.h"
#include "GeomDataLagrange.h"
#include "SolutionData.h"
#include "Functions.h"
#include "QuadratureUtil.h"
#include "myConstants.h"

#include "stabilisationRoutines.h"

#include "KimMoinFlow.h"

using namespace std;

extern ComputerTime       computerTime;
extern MpapTime mpapTime;


LagrangeElem2DStokesQuad4Node::LagrangeElem2DStokesQuad4Node()
{
  degree = 1;
  npElem = 4;
  nlbf   = 4;
  ndof   = 3;
  nsize  = nlbf*ndof;

  if (debug) cout << " constructor LagrangeElem2DStokesQuad4Node\n\n";
}

LagrangeElem2DStokesQuad4Node::~LagrangeElem2DStokesQuad4Node()
{
  if (debug) cout << " destructor LagrangeElem2DStokesQuad4Node\n\n";
}


void LagrangeElem2DStokesQuad4Node::prepareElemData()
{
  LagrangeElement::prepareElemData();

  return;
}


void LagrangeElem2DStokesQuad4Node::prepareElemData2()
{
  return;
}



int LagrangeElem2DStokesQuad4Node::calcLoadVector()
{
  return 0;
}



int LagrangeElem2DStokesQuad4Node::calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal, bool firstIter)
{
    //Stokes2DEx1 analy;
    //Kovasznay  analy;
    //analy.SetPressure(0.0);

    int ii, jj, gp1, gp2, TI, TIp1, TIp2, count, TJ, TJp1, TJp2;

    double  uu, vv, Jac, dvol, b1, b2, b3, b4, b5, b6, b7, b8, xx, yy, acceFact, HH, dist, rho, mu;
    double  pres, Da, Db, af, am, d1, c1, muTaf, rad, urdr, urdr2, tau[3], volume, dt, tCur,  CI=40.0;
    double  fact, fact1, fact2, bb1, bb2, param[2], h2, stabParam, *elmDat, xNode[4], yNode[4];

    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    VectorXd  N(nlbf), dN_dx(nlbf), dN_dy(nlbf);

    VectorXd  res(3), dp(2), Du(2), vel(2), velDot(2), force(2), res2(2), rStab(2);
    MatrixXd  F(2,2), stress(2,2), Dj(2,3);
    Dj.setZero();

    VectorXd  velPrev(3), velTemp(3);

    elmDat = &(SolnData->ElemProp[elmType].data[0]);
    //matDat = &(SolnData->MatlProp[matType].data[0]);

    rho = elmDat[4];
    mu  = elmDat[5];

    KimMoinFlow  analy(rho, mu, 1.0);

    bool axsy = ((int)elmDat[2] == 1);
    //cout << " axsy = " << axsy << endl;

    dt = mpapTime.dt;
    af = SolnData->td(2);
    am = SolnData->td(3);
    acceFact = am*SolnData->td(9);
    muTaf = mu*af;

    tCur = mpapTime.cur - (1.0-af)*dt;

    volume = 0.5*( (xNode[0]-xNode[3])*(yNode[1]-yNode[2]) + (xNode[1]-xNode[2])*(yNode[3]-yNode[0]) );
    //cout << " volume = " << volume << endl;
    h2 = 4.0*volume/PI;

    stabParam = h2/(12.0*mu);
    tau[0] = elmDat[8]*stabParam;      // SUPG
    tau[1] = elmDat[9]*stabParam;//rho;  // PSPG
    tau[2] = elmDat[10]*stabParam*rho; // LSIC

    MatrixXd  matB(2,2), matBinv(2,2), matG(3,3), matJ(2,2), matJinv(2,2);

    double  hx = GeomData->NodePosOrig[nodeNums[1]][0] - GeomData->NodePosOrig[nodeNums[0]][0];
    double  hy = hx;

    matB.setZero();
    matB(0,0) = hx*0.5;
    matB(1,1) = hy*0.5;

    matBinv = matB.inverse();

    matG.setZero();
    matG(0,0) = matBinv(0,0);
    matG(0,1) = matBinv(0,1);
    matG(1,0) = matBinv(1,0);
    matG(1,1) = matBinv(1,1);

    matG = matG.transpose() * matG;

    int nGP1 = 2;
    int nGP2 = 2;

    if(Klocal.rows() != nsize)
    {
      Klocal.resize(nsize, nsize);
      Flocal.resize(nsize);
    }

    Klocal.setZero();
    Flocal.setZero();
    
    //cout << " AAAAAAAAAA " << endl;
    //cout << nGP1 << '\t' << nGP2 << endl;

    vector<double>  gausspoints1, gaussweights1;

    getGaussPoints1D(nGP1, gausspoints1, gaussweights1);


    for(gp2=0;gp2<nGP2;gp2++)
    {
        param[1] = gausspoints1[gp2];

    for(gp1=0;gp1<nGP1;gp1++)
    {
          param[0] = gausspoints1[gp1];

          // basis functions with respect to original configuration
          //GeomData->computeBasisFunctions2D(0, 2, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);
          
          // basis functions with respect to current configuration for ALE formulation
          GeomData->computeBasisFunctions2D(1, 2, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

          dvol = gaussweights1[gp2] * gaussweights1[gp1] * Jac;

          //xx = computeGeomOrig(0, N);
          //yy = computeGeomOrig(1, N);
          xx = yy= 0.0;
          for(ii=0;ii<nlbf;ii++)
          {
            xx += N[ii]*xNode[ii];
            yy += N[ii]*yNode[ii];
          }

          vel(0) = computeValueCur(0, N);
          vel(1) = computeValueCur(1, N);

          velPrev(0) = computeValuePrev(0, N);
          velPrev(1) = computeValuePrev(1, N);

          F(0,0) = computeValueCur(0, dN_dx);
          F(0,1) = computeValueCur(0, dN_dy);
          F(1,0) = computeValueCur(1, dN_dx);
          F(1,1) = computeValueCur(1, dN_dy);
          
          Du.setZero();
          
          pres   = computeValueCur(2, N);
          dp(0)  = computeValueCur(2, dN_dx);
          dp(1)  = computeValueCur(2, dN_dy);

          velDot(0) = computeValueDotCur(0, N);
          velDot(1) = computeValueDotCur(1, N);

          // this is pseudo-stress
          stress = mu*F;
          stress(0,0) -= pres;
          stress(1,1) -= pres;

          //cout << xx << '\t' << yy << endl;
          force.setZero();
          force(0) = analy.computeForce(0, xx, yy, tCur);
          force(1) = analy.computeForce(1, xx, yy, tCur);
          //cout << force(0) << '\t' << force(1) << endl;

          res2(0) = rho*(velDot(0) - force(0)) ;
          res2(1) = rho*(velDot(1) - force(1)) ;

          rStab(0) = res2(0) - mu*Du(0) + dp(0) ;
          rStab(1) = res2(1) - mu*Du(1) + dp(1) ;

          velTemp(0) = velPrev(0);
          velTemp(1) = velPrev(1);
          velTemp(2) = 0.0;

          //evaluateStabParams_algo2(&velTemp(0), h, rho, mu, dt,  beta, tau);

          evaluateStabParams_algo3(velTemp, matG, dt, rho, mu, CI, tau);

          tau[0] *= elmDat[8];  // SUPG
          tau[1] *= elmDat[9];  // PSPG
          tau[2] *= elmDat[10]; // LSIC


          if(axsy)
          {
            //rad = xx;
            //urdr  = vel(0)/rad;

            rad = yy;
            urdr  = vel(1)/rad;

            urdr2 = urdr/rad;
            dvol *= (2.0*PI*rad);
          }

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

            c1 = (rho*acceFact)*b4;

            for(jj=0;jj<nlbf;jj++)
            {
              TJ   = ndof*jj;
              TJp1 = TJ+1;
              TJp2 = TJ+2;

              // time acceleration term
              fact2 = rho*acceFact*N(jj);

              // time acceleration term
              fact = b4*fact2 ;

              // diffusion term
              fact += ( b5*dN_dx(jj)+b6*dN_dy(jj) );

              Klocal(TI,   TJ)   += fact;
              Klocal(TIp1, TJp1) += fact;

              // pressure term
              Klocal(TI,   TJp2) -= (b1*af*N(jj));
              Klocal(TIp1, TJp2) -= (b2*af*N(jj));

              // continuity equation
              Klocal(TIp2, TJ)   -= (b8*dN_dx(jj));
              Klocal(TIp2, TJp1) -= (b8*dN_dy(jj));

              // PSPG stabilisation terms
              //fact2 -= muTaf*d2N(jj);

              Dj(0,0) = fact2;
              Dj(0,1) = 0.0;
              Dj(0,2) = af*dN_dx(jj);
              Dj(1,0) = 0.0;
              Dj(1,1) = fact2;
              Dj(1,2) = af*dN_dy(jj);

              // PSPG
              Klocal(TIp2, TJ)   -= (b1*Dj(0,0) + b2*Dj(1,0))*tau[1];
              Klocal(TIp2, TJp1) -= (b1*Dj(0,1) + b2*Dj(1,1))*tau[1];
              Klocal(TIp2, TJp2) -= (b1*Dj(0,2) + b2*Dj(1,2))*tau[1];

              if(axsy)
              {
                  // diffusion term
                  //Klocal(TI, TJ)     += (mu*b4*af*N(jj)/rad/rad);
                  //Klocal(TI, TJp2)   -= (b4*N(jj)/rad);

                  // continuity equation
                  //Klocal(TIp2, TJ)   -= (b4*af*N(jj)/rad);

                  // diffusion term
                  Klocal(TIp1, TJp1)   += (mu*b4*af*N(jj)/rad/rad);
                  Klocal(TIp1, TJp2)   -= (b4*N(jj)/rad);

                  // continuity equation
                  Klocal(TIp2, TJp1)   -= (b4*af*N(jj)/rad);

              }
            }

            Flocal(TI)   -= (b4*res2(0) + b1*stress(0,0) + b2*stress(0,1) );
            Flocal(TIp1) -= (b4*res2(1) + b1*stress(1,0) + b2*stress(1,1) );
            Flocal(TIp2) += (b4*F.trace());

            // PSPG stabilisation terms
            Flocal(TIp2) += (tau[1]*(b1*rStab(0)+b2*rStab(1)));

            if(axsy)
            {
                //Flocal(TI)   += (-b4*(mu*vel(0)/rad/rad));
                //Flocal(TI)   += (b4*pres/rad);
                //Flocal(TIp2) += (b4*vel(0)/rad);

                Flocal(TIp1)   += (-b4*(mu*vel(1)/rad/rad));
                Flocal(TIp1)   += (b4*pres/rad);
                Flocal(TIp2)   += (b4*vel(1)/rad);
            }
          }
  }//gp1
  }//gp2

  return 0;
}



int LagrangeElem2DStokesQuad4Node::calcInternalForces()
{
  return 0;
}



void LagrangeElem2DStokesQuad4Node::discreteContourplot(int vartype, int varindex, int index, int nCol, double umin, double umax)
{
  return;
}


void LagrangeElem2DStokesQuad4Node::projectToKnots(bool extrapolateFlag, int vartype, int varindex, int index)
{
  return;
}


void LagrangeElem2DStokesQuad4Node::projectStress(int varindex, double* outval)
{
  return;
}



void LagrangeElem2DStokesQuad4Node::projectStrain(int vartype, int varindex, double* outval)
{
  return;
}



void LagrangeElem2DStokesQuad4Node::projectIntVar(int index, double* outval)
{
  return;
}


int LagrangeElem2DStokesQuad4Node::calcOutput(double u1, double v1)
{
  return 0;
}



void LagrangeElem2DStokesQuad4Node::toPostprocess(int vartype, int varindex, int type, SparseMatrixXd&  coeffMat, VectorXd& rhsVec)
{
  return;
}


