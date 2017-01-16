
#include "LagrangeElem2DNavierStokesQuad4Node.h"

#include "Debug.h"
#include "MpapTime.h"
#include "ComputerTime.h"
#include "GeomDataLagrange.h"
#include "SolutionData.h"
#include "Functions.h"
#include "QuadratureUtil.h"

#include "stabilisationRoutines.h"

#include "KimMoinFlow.h"

using namespace std;

extern ComputerTime       computerTime;
extern MpapTime mpapTime;


LagrangeElem2DNavierStokesQuad4Node::LagrangeElem2DNavierStokesQuad4Node()
{
  degree = 1;
  npElem = 4;
  nlbf   = 4;
  ndof   = 3;
  nsize  = nlbf*ndof;

  if (debug) cout << " constructor LagrangeElem2DNavierStokesQuad4Node\n\n";
}

LagrangeElem2DNavierStokesQuad4Node::~LagrangeElem2DNavierStokesQuad4Node()
{
  if (debug) cout << " destructor LagrangeElem2DNavierStokesQuad4Node\n\n";
}


void LagrangeElem2DNavierStokesQuad4Node::prepareElemData()
{
  LagrangeElement::prepareElemData();

  //Klocal.resize(nsize, nsize);
  //Flocal.resize(nsize);

  return;
}


void LagrangeElem2DNavierStokesQuad4Node::prepareElemData2()
{
  return;
}



void LagrangeElem2DNavierStokesQuad4Node::resetMatrixAndVector()
{
  return;
}

int LagrangeElem2DNavierStokesQuad4Node::calcLoadVector()
{
  return 0;
}



/*
int LagrangeElem2DNavierStokesQuad4Node::calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal)
{
  // Fully-implicit formulation
  
  //   char fct[] = "LagrangeElem2DNavierStokesQuad4Node::calcStiffnessAndResidual";
//   computerTime.go(fct);

    //Stokes2DEx1 analy;
    //Kovasznay  analy;
    //analy.SetPressure(0.0);

    int ii, jj, gp1, gp2, TI, TIp1, TIp2, count, TJ, TJp1, TJp2;

    double  uu, vv, Jac, dvol, b1, b2, b3, b4, b5, b6, b7, b8, xx, yy, acceFact, HH, dist, rho, mu;
    double  pres, Da, Db, af, am, d1, c1, muTaf, rad, urdr, urdr2, conv, tau[3], volume, CI=40.0;
    double  fact, fact1, fact2, bb1, bb2, param[2], h2, stabParam, *elmDat, x[4], y[4], dt, tCur;

    for(ii=0;ii<npElem;ii++)
    {
      x[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      y[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    VectorXd  N(nlbf), dN_dx(nlbf), dN_dy(nlbf);
    VectorXd  res(3), res2(2), dp(2), Du(2), vel(2), velDot(2), force(2), gradTvel(2), rStab(3);
    MatrixXd  Dj(2, 3), grad(2,2), gradN(2,2), stress(2,2), gradPrev(2,2);
    Dj.setZero();
    VectorXd  velPrev(3), velTemp(3);


    MatrixXd  matB(2,2), matBinv(2,2), matG(3,3), matJ(2,2), matJinv(2,2);

    elmDat = &(SolnData->ElemProp.data[0]);
    //matDat = &(SolidSolnData->MatlProp.data[0]);

    rho = elmDat[4];
    mu  = elmDat[5];

    KimMoinFlow  analy(rho, mu, 1.0);

    bool axsy = false;
    axsy = ((int)elmDat[2] == 1);

    dt = mpapTime.dt;
    af = SolnData->td(2);
    am = SolnData->td(3);
    acceFact = rho*am*SolnData->td(9);
    muTaf = mu*af;

    tCur = mpapTime.cur - (1.0-af)*dt;

    volume = 0.5*( (x[0]-x[3])*(y[1]-y[2]) + (x[1]-x[2])*(y[3]-y[0]) );
    //cout << " volume = " << volume << endl;
    h2 = 4.0*volume/PI;

    stabParam = h2/(4.0*mu);
    tau[0] = elmDat[8]*stabParam;      // SUPG
    tau[1] = elmDat[9]*stabParam;//rho;  // PSPG
    tau[2] = elmDat[10]*stabParam*rho; // LSIC

    double  hx = GeomData->NodePosOrig[nodeNums[1]][0] - GeomData->NodePosOrig[nodeNums[0]][0];
    double  hy = hx;

    //CI = 1.0/hx;

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


    nGP1 = nGP2 = 2;

    Klocal.setZero();
    Flocal.setZero();
    
    //cout << " AAAAAAAAAA " << endl;
    //cout << nGP1 << '\t' << nGP2 << endl;

    for(gp2=0;gp2<nGP2;gp2++)
    {
        param[1] = GeomData->gausspoints2[gp2];

    for(gp1=0;gp1<nGP1;gp1++)
    {
          param[0] = GeomData->gausspoints1[gp1];

          GeomData->computeBasisFunctions2D(0, 2, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

          dvol = GeomData->gaussweights2[gp2] * GeomData->gaussweights1[gp1] * Jac;

          xx = computeGeomOrig(0, N);
          yy = computeGeomOrig(1, N);

          vel(0) = computeValueCur(0, N);
          vel(1) = computeValueCur(1, N);

          velPrev(0) = computeValuePrev(0, N);
          velPrev(1) = computeValuePrev(1, N);

          grad(0,0) = computeValueCur(0, dN_dx);
          grad(0,1) = computeValueCur(0, dN_dy);
          grad(1,0) = computeValueCur(1, dN_dx);
          grad(1,1) = computeValueCur(1, dN_dy);

          Du.setZero();

          pres   = computeValue(2, N);
          pres   = computeValueCur(2, N);
          dp(0)  = computeValue(2, dN_dx);
          dp(1)  = computeValue(2, dN_dy);

          velDot(0) = computeValueDotCur(0, N);
          velDot(1) = computeValueDotCur(1, N);

          // this is pseudo-stress
          stress = mu*grad;
          stress(0,0) -= pres;
          stress(1,1) -= pres;

          //cout << xx << '\t' << yy << endl;

          force.setZero();
          force(0) = analy.computeForce(0, xx, yy, tCur);
          force(1) = analy.computeForce(1, xx, yy, tCur);
          //cout << force(0) << '\t' << force(1) << endl;

          if(axsy)
          {
            rad = xx;

            urdr  = vel(0)/rad;
            urdr2 = urdr/rad;
            dvol *= (2.0*PI*rad);
          }

          gradTvel = grad*vel;

          res2(0) = rho*(velDot(0) + gradTvel(0) - force(0)) ;
          res2(1) = rho*(velDot(1) + gradTvel(1) - force(1)) ;

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
             
             Da = (vel(0)*b1 + vel(1)*b2)*tau[0];

             for(jj=0;jj<nlbf;jj++)
             {
               TJ   = ndof*jj;
               TJp1 = TJ+1;
               TJp2 = TJ+2;

               fact2 = rho*acceFact*N(jj);
               
               // time acceleration term
               fact = b4*fact2 ;

               // diffusion term
               fact += b5*dN_dx(jj)+b6*dN_dy(jj);

               Klocal(TI,   TJ)   += fact;
               Klocal(TIp1, TJp1) += fact;

               // convection term

               gradN = grad*(rho*N(jj));

               Db = rho*(vel(0)*dN_dx(jj) + vel(1)*dN_dy(jj));
               
               gradN(0,0) += Db;
               gradN(1,1) += Db;

               Klocal(TI,   TJ)   += (b8*gradN(0,0));
               Klocal(TI,   TJp1) += (b8*gradN(0,1));
               Klocal(TIp1, TJ)   += (b8*gradN(1,0));
               Klocal(TIp1, TJp1) += (b8*gradN(1,1));

               // pressure term
               Klocal(TI,   TJp2) -= (b1*af*N(jj));
               Klocal(TIp1, TJp2) -= (b2*af*N(jj));

               // continuity equation
               Klocal(TIp2, TJ)   += (b8*dN_dx(jj));
               Klocal(TIp2, TJp1) += (b8*dN_dy(jj));
               Klocal(TIp2, TJp2) += 0.0; //

               // SUPG and PSPG stabilisation terms
               //fact2 -= mu*d2N(jj);
               fact2 -= 0.0;
               
               gradN *= af;
               
               Dj(0,0) = gradN(0,0) + fact2;
               Dj(0,1) = gradN(0,1);
               Dj(0,2) = dN_dx(jj);
               Dj(1,0) = gradN(1,0);
               Dj(1,1) = gradN(1,1) + fact2;
               Dj(1,2) = dN_dy(jj);

               // SUPG
               Klocal(TI, TJ)     += Da*Dj(0,0);
               Klocal(TI, TJp1)   += Da*Dj(0,1);
               Klocal(TI, TJp2)   += Da*Dj(0,2);

               Klocal(TIp1, TJ)   += Da*Dj(1,0);
               Klocal(TIp1, TJp1) += Da*Dj(1,1);
               Klocal(TIp1, TJp2) += Da*Dj(1,2);

               Klocal(TI,   TJ)   += ( (tau[0]*af) * b1 * rStab(0) * N(jj) );
               Klocal(TI,   TJp1) += ( (tau[0]*af) * b2 * rStab(0) * N(jj) );
               Klocal(TIp1, TJ)   += ( (tau[0]*af) * b1 * rStab(1) * N(jj) );
               Klocal(TIp1, TJp1) += ( (tau[0]*af) * b2 * rStab(1) * N(jj) );

               // PSPG
               Klocal(TIp2, TJ)   += (b1*Dj(0,0) + b2*Dj(1,0))*tau[1];
               Klocal(TIp2, TJp1) += (b1*Dj(0,1) + b2*Dj(1,1))*tau[1];
               Klocal(TIp2, TJp2) += (b1*Dj(0,2) + b2*Dj(1,2))*tau[1];
               
               // LSIC stabilisation

               fact2 = rho*af*tau[2];
               Klocal(TI,   TJ)   += (b1*dN_dx(jj))*fact2;
               Klocal(TI,   TJp1) += (b1*dN_dy(jj))*fact2;

               Klocal(TIp1, TJ)   += (b2*dN_dx(jj))*fact2;
               Klocal(TIp1, TJp1) += (b2*dN_dy(jj))*fact2;

               if(axsy)
               {
                  // diffusion term
                  Klocal(TI, TJ)     += (mu*b4*N(jj)/rad/rad);
                  Klocal(TI, TJp2)   -= (b4*N(jj)/rad);

                  // continuity equation
                  Klocal(TIp2, TJ)   -= (b4*N(jj)/rad);
               }
             }

             Flocal(TI)   -= (b4*res2(0) + b1*stress(0,0) + b2*stress(0,1) );
             Flocal(TIp1) -= (b4*res2(1) + b1*stress(1,0) + b2*stress(1,1) );
             Flocal(TIp2) -= (b4*grad.trace());

             // SUPG stabilisation terms
             Flocal(TI)   -= Da*rStab(0);
             Flocal(TIp1) -= Da*rStab(1);
             
             // PSPG stabilisation terms
             Flocal(TIp2) -= (tau[1]*(b1*rStab(0)+b2*rStab(1)));

             // LSIC stabilisation terms
             fact2 = tau[2]*rho*grad.trace();

             Flocal(TI)   -= b1*fact2;
             Flocal(TIp1) -= b2*fact2;

             if(axsy)
             {
                Flocal(TI)   += (-b4*(mu*vel(0)/rad/rad));
                Flocal(TI)   += (b4*pres/rad);
                Flocal(TIp2) += (b4*vel(0)/rad);
             }
          }
  }//gp1
  }//gp2

  return 0;
}
*/


//
int LagrangeElem2DNavierStokesQuad4Node::calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal)
{
  // Semi-implicit formulation - type A
  
  //   char fct[] = "LagrangeElem2DNavierStokesQuad4Node::calcStiffnessAndResidual";
//   computerTime.go(fct);

    //Stokes2DEx1 analy;
    //Kovasznay  analy;
    //analy.SetPressure(0.0);

    int ii, jj, gp1, gp2, TI, TIp1, TIp2, count, TJ, TJp1, TJp2;

    double  uu, vv, Jac, dvol, b1, b2, b3, b4, b5, b6, b7, b8, xx, yy, acceFact, HH, dist, rho, mu;
    double  pres, Da, Db, af, am, d1, c1, muTaf, rad, urdr, urdr2, conv, tau[3], volume, CI=40.0;
    double  fact, fact1, fact2, bb1, bb2, param[2], h2, stabParam, *elmDat, x[4], y[4], dt, tCur;

    for(ii=0;ii<npElem;ii++)
    {
      x[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      y[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    VectorXd  N(nlbf), dN_dx(nlbf), dN_dy(nlbf);
    VectorXd  res(3), res2(2), dp(2), Du(2), vel(2), velDot(2), force(2), gradTvel(2), rStab(3);
    MatrixXd  Dj(2, 3), grad(2,2), gradN(2,2), stress(2,2), gradPrev(2,2);
    Dj.setZero();
    VectorXd  velPrev(2), velExt(2), velTemp(3);


    MatrixXd  matB(2,2), matBinv(2,2), matG(3,3), matJ(2,2), matJinv(2,2);

    elmDat = &(SolnData->ElemProp[elmType].data[0]);
    //matDat = &(SolidSolnData->MatlProp[matType].data[0]);

    rho = elmDat[4];
    mu  = elmDat[5];

    KimMoinFlow  analy(rho, mu, 1.0);

    bool axsy = false;
    axsy = ((int)elmDat[2] == 1);

    dt = mpapTime.dt;
    af = SolnData->td(2);
    am = SolnData->td(3);
    acceFact = rho*am*SolnData->td(9);
    muTaf = mu*af;

    tCur = mpapTime.cur - (1.0-af)*dt;

    volume = 0.5*( (x[0]-x[3])*(y[1]-y[2]) + (x[1]-x[2])*(y[3]-y[0]) );
    //cout << " volume = " << volume << endl;
    h2 = 4.0*volume/PI;

    stabParam = h2/(4.0*mu);
    tau[0] = elmDat[8]*stabParam;      // SUPG
    tau[1] = elmDat[9]*stabParam;//rho;  // PSPG
    tau[2] = elmDat[10]*stabParam*rho; // LSIC

    double  hx = GeomData->NodePosOrig[nodeNums[1]][0] - GeomData->NodePosOrig[nodeNums[0]][0];
    double  hy = hx;

    //CI = 1.0/hx;

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

    for(gp2=0;gp2<nGP2;gp2++)
    {
        param[1] = GeomData->gausspoints2[gp2];

    for(gp1=0;gp1<nGP1;gp1++)
    {
          param[0] = GeomData->gausspoints1[gp1];

          GeomData->computeBasisFunctions2D(0, 2, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

          dvol = GeomData->gaussweights2[gp2] * GeomData->gaussweights1[gp1] * Jac;

          xx = computeGeomOrig(0, N);
          yy = computeGeomOrig(1, N);

          vel(0) = computeValueCur(0, N);
          vel(1) = computeValueCur(1, N);

          velPrev(0) = computeValuePrev(0, N);
          velPrev(1) = computeValuePrev(1, N);

          velExt(0) = computeValueExtrap(0, N);
          velExt(1) = computeValueExtrap(1, N);

          grad(0,0) = computeValueCur(0, dN_dx);
          grad(0,1) = computeValueCur(0, dN_dy);
          grad(1,0) = computeValueCur(1, dN_dx);
          grad(1,1) = computeValueCur(1, dN_dy);

          gradPrev(0,0) = computeValuePrev(0, dN_dx);
          gradPrev(0,1) = computeValuePrev(0, dN_dy);
          gradPrev(1,0) = computeValuePrev(1, dN_dx);
          gradPrev(1,1) = computeValuePrev(1, dN_dy);

          Du.setZero();

          pres   = computeValueCur(2, N);
          dp(0)  = computeValueCur(2, dN_dx);
          dp(1)  = computeValueCur(2, dN_dy);

          velDot(0) = computeValueDotCur(0, N);
          velDot(1) = computeValueDotCur(1, N);

          // this is pseudo-stress
          stress = mu*grad;
          stress(0,0) -= pres;
          stress(1,1) -= pres;

          //cout << xx << '\t' << yy << endl;

          force.setZero();
          force(0) = analy.computeForce(0, xx, yy, tCur);
          force(1) = analy.computeForce(1, xx, yy, tCur);
          //cout << force(0) << '\t' << force(1) << endl;

          if(axsy)
          {
            rad = xx;

            urdr  = vel(0)/rad;
            urdr2 = urdr/rad;
            dvol *= (2.0*PI*rad);
          }

          velExt = af*velExt + (1.0-af)*velPrev;

          gradTvel = grad*velExt;

          res2(0) = rho*(velDot(0) + gradTvel(0) - force(0)) ;
          res2(1) = rho*(velDot(1) + gradTvel(1) - force(1)) ;

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
             
             //Da = (velPrev(0)*b1 + velPrev(1)*b2)*tau[0];
             Da = (velExt(0)*b1 + velExt(1)*b2)*tau[0];

             for(jj=0;jj<nlbf;jj++)
             {
               TJ   = ndof*jj;
               TJp1 = TJ+1;
               TJp2 = TJ+2;

               fact2 = rho*acceFact*N(jj);
               
               // time acceleration term
               fact = b4*fact2 ;

               // diffusion term
               fact += b5*dN_dx(jj)+b6*dN_dy(jj);

               Klocal(TI,   TJ)   += fact;
               Klocal(TIp1, TJp1) += fact;

               // convection term

               //Db = rho*(velPrev(0)*dN_dx(jj) + velPrev(1)*dN_dy(jj));
               Db = rho*(velExt(0)*dN_dx(jj) + velExt(1)*dN_dy(jj));

               Klocal(TI,   TJ)   += (b8*Db);
               Klocal(TIp1, TJp1) += (b8*Db);

               // pressure term
               Klocal(TI,   TJp2) -= (af*b1*N(jj));
               Klocal(TIp1, TJp2) -= (af*b2*N(jj));

               // continuity equation
               Klocal(TIp2, TJ)   += (b8*dN_dx(jj));
               Klocal(TIp2, TJp1) += (b8*dN_dy(jj));

               // SUPG and PSPG stabilisation terms
               //fact2 -= mu*d2N(jj);
               fact2 -= 0.0;

               Db *= af;

               Dj(0,0) = Db + fact2;
               Dj(0,1) = 0.0;
               Dj(0,2) = af*dN_dx(jj);
               Dj(1,0) = 0.0;
               Dj(1,1) = Db + fact2;
               Dj(1,2) = af*dN_dy(jj);

               // SUPG
               Klocal(TI, TJ)     += Da*Dj(0,0);
               Klocal(TI, TJp1)   += Da*Dj(0,1);
               Klocal(TI, TJp2)   += Da*Dj(0,2);

               Klocal(TIp1, TJ)   += Da*Dj(1,0);
               Klocal(TIp1, TJp1) += Da*Dj(1,1);
               Klocal(TIp1, TJp2) += Da*Dj(1,2);

               // PSPG
               Klocal(TIp2, TJ)   += (b1*Dj(0,0) + b2*Dj(1,0))*tau[1];
               Klocal(TIp2, TJp1) += (b1*Dj(0,1) + b2*Dj(1,1))*tau[1];
               Klocal(TIp2, TJp2) += (b1*Dj(0,2) + b2*Dj(1,2))*tau[1];
               
               // LSIC stabilisation

               fact2 = rho*af*tau[2];
               Klocal(TI,   TJ)   += (b1*dN_dx(jj))*fact2;
               Klocal(TI,   TJp1) += (b1*dN_dy(jj))*fact2;

               Klocal(TIp1, TJ)   += (b2*dN_dx(jj))*fact2;
               Klocal(TIp1, TJp1) += (b2*dN_dy(jj))*fact2;

               if(axsy)
               {
                  // diffusion term
                  Klocal(TI, TJ)     += (mu*b4*N(jj)/rad/rad);
                  Klocal(TI, TJp2)   -= (b4*N(jj)/rad);

                  // continuity equation
                  Klocal(TIp2, TJ)   -= (b4*N(jj)/rad);
               }
             }

             Flocal(TI)   -= (b4*res2(0) + b1*stress(0,0) + b2*stress(0,1) );
             Flocal(TIp1) -= (b4*res2(1) + b1*stress(1,0) + b2*stress(1,1) );
             Flocal(TIp2) -= (b4*grad.trace());

             // SUPG stabilisation terms
             Flocal(TI)   -= Da*rStab(0);
             Flocal(TIp1) -= Da*rStab(1);
             
             // PSPG stabilisation terms
             Flocal(TIp2) -= (tau[1]*(b1*rStab(0)+b2*rStab(1)));

             // LSIC stabilisation terms

             fact2 = tau[2]*rho*grad.trace();

             Flocal(TI)   -= b1*fact2;
             Flocal(TIp1) -= b2*fact2;

             if(axsy)
             {
                Flocal(TI)   += (-b4*(mu*vel(0)/rad/rad));
                Flocal(TI)   += (b4*pres/rad);
                Flocal(TIp2) += (b4*vel(0)/rad);
             }
          }
  }//gp1
  }//gp2

  return 0;
}
//




/*
int LagrangeElem2DNavierStokesQuad4Node::calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal)
{
  // Semi-implicit formulation - type B
  
  //   char fct[] = "LagrangeElem2DNavierStokesQuad4Node::calcStiffnessAndResidual";
//   computerTime.go(fct);

    //Stokes2DEx1 analy;
    //Kovasznay  analy;
    //analy.SetPressure(0.0);

    int ii, jj, gp1, gp2, TI, TIp1, TIp2, count, TJ, TJp1, TJp2;

    double  uu, vv, Jac, dvol, b1, b2, b3, b4, b5, b6, b7, b8, xx, yy, acceFact, HH, dist, rho, mu;
    double  pres, Da, Db, af, am, d1, c1, muTaf, rad, urdr, urdr2, conv, tau[3], volume, CI=40.0;
    double  fact, fact1, fact2, bb1, bb2, param[2], h2, stabParam, *elmDat, x[4], y[4], dt, tCur;

    for(ii=0;ii<npElem;ii++)
    {
      x[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      y[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    VectorXd  N(nlbf), dN_dx(nlbf), dN_dy(nlbf);
    VectorXd  res(3), res2(2), dp(2), Du(2), vel(2), velDot(2), force(2), gradTvel(2), rStab(3);
    MatrixXd  Dj(2, 3), grad(2,2), gradN(2,2), stress(2,2), gradPrev(2,2);
    Dj.setZero();
    VectorXd  velPrev(2), velTemp(3);


    MatrixXd  matB(2,2), matBinv(2,2), matG(3,3), matJ(2,2), matJinv(2,2);

    elmDat = &(SolnData->ElemProp.data[0]);
    //matDat = &(SolidSolnData->MatlProp.data[0]);

    rho = elmDat[4];
    mu  = elmDat[5];

    KimMoinFlow  analy(rho, mu, 1.0);

    bool axsy = false;
    axsy = ((int)elmDat[2] == 1);

    dt = mpapTime.dt;
    af = SolnData->td(2);
    am = SolnData->td(3);
    acceFact = rho*am*SolnData->td(9);
    muTaf = mu*af;

    tCur = mpapTime.cur - (1.0-af)*dt;

    volume = 0.5*( (x[0]-x[3])*(y[1]-y[2]) + (x[1]-x[2])*(y[3]-y[0]) );
    //cout << " volume = " << volume << endl;
    h2 = 4.0*volume/PI;

    stabParam = h2/(4.0*mu);
    tau[0] = elmDat[8]*stabParam;      // SUPG
    tau[1] = elmDat[9]*stabParam;//rho;  // PSPG
    tau[2] = elmDat[10]*stabParam*rho; // LSIC

    double  hx = GeomData->NodePosOrig[nodeNums[1]][0] - GeomData->NodePosOrig[nodeNums[0]][0];
    double  hy = hx;

    //CI = 1.0/hx;

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


    nGP1 = nGP2 = 2;

    Klocal.setZero();
    Flocal.setZero();
    
    //cout << " AAAAAAAAAA " << endl;
    //cout << nGP1 << '\t' << nGP2 << endl;

    for(gp2=0;gp2<nGP2;gp2++)
    {
        param[1] = GeomData->gausspoints2[gp2];

    for(gp1=0;gp1<nGP1;gp1++)
    {
          param[0] = GeomData->gausspoints1[gp1];

          GeomData->computeBasisFunctions2D(0, 2, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

          dvol = GeomData->gaussweights2[gp2] * GeomData->gaussweights1[gp1] * Jac;

          xx = computeGeomOrig(0, N);
          yy = computeGeomOrig(1, N);

          vel(0) = computeValueCur(0, N);
          vel(1) = computeValueCur(1, N);

          velPrev(0) = computeValuePrev(0, N);
          velPrev(1) = computeValuePrev(1, N);

          grad(0,0) = computeValueCur(0, dN_dx);
          grad(0,1) = computeValueCur(0, dN_dy);
          grad(1,0) = computeValueCur(1, dN_dx);
          grad(1,1) = computeValueCur(1, dN_dy);

          gradPrev(0,0) = computeValuePrev(0, dN_dx);
          gradPrev(0,1) = computeValuePrev(0, dN_dy);
          gradPrev(1,0) = computeValuePrev(1, dN_dx);
          gradPrev(1,1) = computeValuePrev(1, dN_dy);

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
          force(0) = analy.computeForce(0, xx, yy, tCur);
          force(1) = analy.computeForce(1, xx, yy, tCur);
          //cout << force(0) << '\t' << force(1) << endl;

          if(axsy)
          {
            rad = xx;

            urdr  = vel(0)/rad;
            urdr2 = urdr/rad;
            dvol *= (2.0*PI*rad);
          }

          gradTvel = gradPrev*vel + grad*velPrev - gradPrev*velPrev;

          res2(0) = rho*(velDot(0) + gradTvel(0) - force(0)) ;
          res2(1) = rho*(velDot(1) + gradTvel(1) - force(1)) ;

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
             
             Da = rho*(velPrev(0)*b1 + velPrev(1)*b2)*tau[0];

             for(jj=0;jj<nlbf;jj++)
             {
               TJ   = ndof*jj;
               TJp1 = TJ+1;
               TJp2 = TJ+2;

               fact2 = rho*acceFact*N(jj);
               
               // time acceleration term
               fact = b4*fact2 ;

               // diffusion term
               fact += b5*dN_dx(jj)+b6*dN_dy(jj);

               Klocal(TI,   TJ)   += fact;
               Klocal(TIp1, TJp1) += fact;

              // convection term - semi-implicit type B

              gradN = gradPrev*(rho*N(jj));

              Db = rho*(velPrev(0)*dN_dx(jj) + velPrev(1)*dN_dy(jj));

              gradN(0,0) += Db;
              gradN(1,1) += Db;

              Klocal(TI,   TJ)   += (b8*gradN(0,0));
              Klocal(TI,   TJp1) += (b8*gradN(0,1));
              Klocal(TIp1, TJ)   += (b8*gradN(1,0));
              Klocal(TIp1, TJp1) += (b8*gradN(1,1));

               // pressure term
               Klocal(TI,   TJp2) -= (b1*N(jj));
               Klocal(TIp1, TJp2) -= (b2*N(jj));

               // continuity equation
               Klocal(TIp2, TJ)   += (b8*dN_dx(jj));
               Klocal(TIp2, TJp1) += (b8*dN_dy(jj));

               // SUPG and PSPG stabilisation terms
               //fact2 -= mu*d2N(jj);
               fact2 -= 0.0;
               
               gradN *= af;
               
               Dj(0,0) = gradN(0,0) + fact2;
               Dj(0,1) = gradN(0,1);
               Dj(0,2) = dN_dx(jj);
               Dj(1,0) = gradN(1,0);
               Dj(1,1) = gradN(1,1) + fact2;
               Dj(1,2) = dN_dy(jj);

               // SUPG
               Klocal(TI, TJ)     += Da*Dj(0,0);
               Klocal(TI, TJp1)   += Da*Dj(0,1);
               Klocal(TI, TJp2)   += Da*Dj(0,2);

               Klocal(TIp1, TJ)   += Da*Dj(1,0);
               Klocal(TIp1, TJp1) += Da*Dj(1,1);
               Klocal(TIp1, TJp2) += Da*Dj(1,2);

               // PSPG
               Klocal(TIp2, TJ)   += (b1*Dj(0,0) + b2*Dj(1,0))*tau[1];
               Klocal(TIp2, TJp1) += (b1*Dj(0,1) + b2*Dj(1,1))*tau[1];
               Klocal(TIp2, TJp2) += (b1*Dj(0,2) + b2*Dj(1,2))*tau[1];
               
               // LSIC stabilisation

               fact2 = rho*af*tau[2];
               Klocal(TI,   TJ)   += (b1*dN_dx(jj))*fact2;
               Klocal(TI,   TJp1) += (b1*dN_dy(jj))*fact2;

               Klocal(TIp1, TJ)   += (b2*dN_dx(jj))*fact2;
               Klocal(TIp1, TJp1) += (b2*dN_dy(jj))*fact2;

               if(axsy)
               {
                  // diffusion term
                  Klocal(TI, TJ)     += (mu*b4*N(jj)/rad/rad);
                  Klocal(TI, TJp2)   -= (b4*N(jj)/rad);

                  // continuity equation
                  Klocal(TIp2, TJ)   -= (b4*N(jj)/rad);
               }
             }

             Flocal(TI)   -= (b4*res2(0) + b1*stress(0,0) + b2*stress(0,1) );
             Flocal(TIp1) -= (b4*res2(1) + b1*stress(1,0) + b2*stress(1,1) );
             Flocal(TIp2) -= (b4*grad.trace());

             // SUPG stabilisation terms
             Flocal(TI)   -= Da*rStab(0);
             Flocal(TIp1) -= Da*rStab(1);
             
             // PSPG stabilisation terms
             Flocal(TIp2) -= (tau[1]*(b1*rStab(0)+b2*rStab(1)));

             // LSIC stabilisation terms
             fact2 = tau[2]*rho*grad.trace();

             Flocal(TI)   -= b1*fact2;
             Flocal(TIp1) -= b2*fact2;

             if(axsy)
             {
                Flocal(TI)   += (-b4*(mu*vel(0)/rad/rad));
                Flocal(TI)   += (b4*pres/rad);
                Flocal(TIp2) += (b4*vel(0)/rad);
             }
          }
  }//gp1
  }//gp2

  return 0;
}
*/




int LagrangeElem2DNavierStokesQuad4Node::calcInternalForces()
{
  return 0;
}



void LagrangeElem2DNavierStokesQuad4Node::discreteContourplot(int vartype, int varindex, int index, int nCol, double umin, double umax)
{
  return;
}


void LagrangeElem2DNavierStokesQuad4Node::projectToKnots(bool extrapolateFlag, int vartype, int varindex, int index)
{
  return;
}


void LagrangeElem2DNavierStokesQuad4Node::projectStress(int varindex, double* outval)
{
  return;
}



void LagrangeElem2DNavierStokesQuad4Node::projectStrain(int vartype, int varindex, double* outval)
{
  return;
}



void LagrangeElem2DNavierStokesQuad4Node::projectIntVar(int index, double* outval)
{
  return;
}


int LagrangeElem2DNavierStokesQuad4Node::calcOutput(double u1, double v1)
{
  return 0;
}



void LagrangeElem2DNavierStokesQuad4Node::toPostprocess(int vartype, int varindex, int type, SparseMatrixXd&  coeffMat, VectorXd& rhsVec)
{
  return;
}


