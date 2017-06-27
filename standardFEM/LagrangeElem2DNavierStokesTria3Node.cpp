
#include "LagrangeElem2DNavierStokesTria3Node.h"

#include "Debug.h"
#include "MpapTime.h"
#include "ComputerTime.h"
#include "GeomDataLagrange.h"
#include "SolutionData.h"
#include "Functions.h"
#include "QuadratureUtil.h"

#include "stabilisationRoutines.h"

#include "BasisFunctionsLagrange.h"

using namespace std;

extern ComputerTime       computerTime;
extern MpapTime mpapTime;


LagrangeElem2DNavierStokesTria3Node::LagrangeElem2DNavierStokesTria3Node()
{
  degree = 1;
  npElem = 3;
  nlbf   = 3;
  ndof   = 3;
  nsize  = nlbf*ndof;

  if (debug) cout << " constructor LagrangeElem2DNavierStokesTria3Node\n\n";
}

LagrangeElem2DNavierStokesTria3Node::~LagrangeElem2DNavierStokesTria3Node()
{
  if (debug) cout << " destructor LagrangeElem2DNavierStokesTria3Node\n\n";
}


void LagrangeElem2DNavierStokesTria3Node::prepareElemData()
{
  LagrangeElement::prepareElemData();

  return;
}


void LagrangeElem2DNavierStokesTria3Node::prepareElemData2()
{
  return;
}

int LagrangeElem2DNavierStokesTria3Node::calcLoadVector()
{
  return 0;
}


//
int LagrangeElem2DNavierStokesTria3Node::calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal)
{
//   char fct[] = "LagrangeElem2DNavierStokes4Node::calcStiffnessAndResidual";
//   computerTime.go(fct);

  // stabilised formulation
  // fully-implicit scheme 

    int ii, jj, gp1, gp2, TI, TIp1, TIp2, count, TJ, TJp1, TJp2;

    double  uu, vv, Jac, dvol, b1, b2, b3, b4, b5, b6, b7, b8, xx, yy, acceFact, rho, mu, CI=1.0;
    double  pres, Da, Db, af, am, d1, c1, muTaf, rad, urdr, urdr2, tau[3], volume;
    double  fact, fact1, fact2, bb1, bb2, param[2], dt, h, h2, stabParam, *elmDat, x[3], y[3], beta[6];
    double  dN_du[2][3];

    beta[0] = 1.0;     beta[1] = 1.0/3.0;
    beta[2] = 30.0;    beta[3] = 0.1;
    beta[4] = 1.0;     beta[5] = 1.0;

    for(ii=0;ii<npElem;ii++)
    {
      x[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      y[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    VectorXd  N(nlbf), dN_dx(nlbf), dN_dy(nlbf), velTemp(3);
    VectorXd  res(3), dp(2), Du(2), vel(2), velDot(2), force(2), res2(2), rStab(2), gradTvel(2), velPrev(2);
    MatrixXd  Dj(2, 3), grad(2,2), gradN(2,2), stress(2,2), matG(3,3), matJ(2,2), matJinv(2,2);

    elmDat = &(SolnData->ElemProp[elmType].data[0]);
    //matDat = &(SolidSolnData->MatlProp.data[0]);

    rho = elmDat[4];
    mu  = elmDat[5];

    bool axsy = false;
    axsy = ((int)elmDat[2] == 1);

    af = SolnData->td(2);
    am = SolnData->td(3);
    acceFact = am*SolnData->td(9);
    muTaf = mu*af;
    dt = mpapTime.dt;

    volume = 0.5*(x[0]*(y[1]-y[2]) + x[1]*(y[2]-y[0]) + x[2]*(y[0]-y[1]));

    //cout << " volume = " << volume << endl;
    h2 = 4.0*volume/PI;
    h = sqrt(h2);

    stabParam = h2/(4.0*mu);
    tau[0] = elmDat[8]*stabParam;  // SUPG
    tau[1] = elmDat[9]*stabParam;  // PSPG
    tau[2] = elmDat[10]*stabParam; // LSIC

    nGP = elmDat[0];

    getGaussPointsTriangle(nGP, GeomData->gausspoints1, GeomData->gausspoints2, GeomData->gaussweights1);

    if(Klocal.rows() != nsize)
    {
      Klocal.resize(nsize, nsize);
      Flocal.resize(nsize);
    }

    Klocal.setZero();
    Flocal.setZero();
    
    //cout << " AAAAAAAAAA " << endl;
    //cout << nGP1 << '\t' << nGP2 << endl;

    for(gp1=0;gp1<nGP;gp1++)
    {
          param[0] = GeomData->gausspoints1[gp1];
          param[1] = GeomData->gausspoints2[gp1];

          GeomData->computeBasisFunctions2D(0, 1, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

          dvol = GeomData->gaussweights1[gp1] * Jac;
          
          LagrangeBasisFunsTria(degree, param[0], param[1], &N(0), dN_du[0], dN_du[1]);
          
          //matJ.setZero();
          //for(ii=0;ii<npElem;ii++)
          //{
            //matJ(0,0) +=  (x[ii] * dN_du[0][ii]) ;
            //matJ(1,0) +=  (x[ii] * dN_du[1][ii]) ;
            //matJ(0,1) +=  (y[ii] * dN_du[0][ii]) ;
            //matJ(1,1) +=  (y[ii] * dN_du[1][ii]) ;
          //}

          //matJinv = matJ.inverse();

          //matG.setZero();
          //matG(0,0) = matJinv(0,0);
          //matG(0,1) = matJinv(0,1);
          //matG(1,0) = matJinv(1,0);
          //matG(1,1) = matJinv(1,1);
          
          //printf("%12.6f \t %12.6f \t %12.6f \n", volume, h, Jac);
          //printf("%12.6f \t %12.6f \t %12.6f \t %12.6f \n", matB(0,0), matB(0,1), matB(1,0), matB(1,1));
          //printf("%12.6f \t %12.6f \t %12.6f \t %12.6f \n", matG(0,0), matG(0,1), matG(1,0), matG(1,1));

          matG = matG.transpose()*matG;

          vel(0) = computeValueCur(0, N);
          vel(1) = computeValueCur(1, N);

          grad(0,0) = computeValueCur(0, dN_dx);
          grad(0,1) = computeValueCur(0, dN_dy);
          grad(1,0) = computeValueCur(1, dN_dx);
          grad(1,1) = computeValueCur(1, dN_dy);

          velPrev(0) = computeValuePrev(0, N);
          velPrev(1) = computeValuePrev(1, N);

          Du.setZero();

          pres   = computeValue(2, N);
          dp(0)  = computeValue(2, dN_dx);
          dp(1)  = computeValue(2, dN_dy);

          velDot(0) = computeValueDotCur(0, N);
          velDot(1) = computeValueDotCur(1, N);
          //velDot.setZero();

          // this is pseudo-stress
          //stress = mu*grad;
          stress = mu*(grad+grad.transpose());
          stress(0,0) -= pres;
          stress(1,1) -= pres;

          //cout << xx << '\t' << yy << endl;
          force.setZero();
          //force(0) = analy.computeForce(0, xx, yy);
          //force(1) = analy.computeForce(1, xx, yy);
          //cout << force(0) << '\t' << force(1) << endl;

          //force(0) = computeForce(0, N);
          //force(1) = computeForce(1, N);

          gradTvel = grad*vel ;

          res2(0) = rho*(velDot(0) + gradTvel(0) - force(0)) ;
          res2(1) = rho*(velDot(1) + gradTvel(1) - force(1)) ;

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

          velTemp(0) = velPrev(0);
          velTemp(1) = velPrev(1);
          velTemp(2) = 0.0;

          //evaluateStabParams_algo1(&velTemp(0), h, rho, mu, dt,  beta, tau);

          //evaluateStabParams_algo2(&velTemp(0), h, rho, mu, dt,  beta, tau);

          //evaluateStabParams_algo3(velTemp, matG, dt, rho, mu, CI, tau);

          //tau[0] *= elmDat[8];  // SUPG
          //tau[1] *= elmDat[9];  // PSPG
          //tau[2] *= elmDat[10]; // LSIC

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

            Da = rho*(vel(0)*b1 + vel(1)*b2)*tau[0];

            for(jj=0;jj<nlbf;jj++)
            {
              TJ   = ndof*jj;
              TJp1 = TJ+1;
              TJp2 = TJ+2;

              fact2 = rho*acceFact*N(jj);

              // time acceleration term
              fact = b4*fact2 ;

              //// diffusion term
              fact += ( b5*dN_dx(jj)+b6*dN_dy(jj) );

              Klocal(TI,   TJ)   += fact;
              Klocal(TIp1, TJp1) += fact;

              Klocal(TI,   TJ)   += ( b5*dN_dx(jj) );
              Klocal(TI,   TJp1) += ( b5*dN_dy(jj) );
              Klocal(TIp1, TJ)   += ( b6*dN_dx(jj) );
              Klocal(TIp1, TJp1) += ( b6*dN_dy(jj) );

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
              Klocal(TI,   TJp2) -= (b1*N(jj));
              Klocal(TIp1, TJp2) -= (b2*N(jj));

              // continuity equation
              Klocal(TIp2, TJ)   += (b8*dN_dx(jj));
              Klocal(TIp2, TJp1) += (b8*dN_dy(jj));

              // SUPG and PSPG stabilisation terms

              gradN *= af;

              Dj(0,0) = gradN(0,0) + fact2;
              Dj(0,1) = gradN(0,1);
              Dj(0,2) = dN_dx(jj);
              Dj(1,0) = gradN(1,0);
              Dj(1,1) = gradN(1,1) + fact2;
              Dj(1,2) = dN_dy(jj);

              if(axsy)
              {
                Dj(0,0) -= muTaf*(dN_dx(jj)/rad - N(jj)/rad/rad);
                Dj(1,1) -= muTaf*(dN_dx(jj)/rad);
              }

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
//



/*
int LagrangeElem2DNavierStokesTria3Node::calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal)
{
//   char fct[] = "LagrangeElem2DNavierStokes4Node::calcStiffnessAndResidual";
//   computerTime.go(fct);

  // stabilised formulation
  // semi-implicit scheme - type A

    int ii, jj, gp1, gp2, TI, TIp1, TIp2, count, TJ, TJp1, TJp2;

    double  uu, vv, Jac, dvol, b1, b2, b3, b4, b5, b6, b7, b8, xx, yy, acceFact, rho, mu, CI=1.0;
    double  pres, Da, Db, af, am, d1, c1, muTaf, rad, urdr, urdr2, tau[3], volume;
    double  fact, fact1, fact2, bb1, bb2, param[2], dt, h, h2, stabParam, *elmDat, beta[6];
    double  dN_du[2][3], xNode[3], yNode[3];

    beta[0] = 1.0;     beta[1] = 1.0/3.0;
    beta[2] = 30.0;    beta[3] = 0.1;
    beta[4] = 1.0;     beta[5] = 1.0;

    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[SolnData->node_map_new_to_old[nodeNums[ii]]][0];
      yNode[ii] = GeomData->NodePosOrig[SolnData->node_map_new_to_old[nodeNums[ii]]][1];
    }


    VectorXd  N(nlbf), dN_dx(nlbf), dN_dy(nlbf), velTemp(3), velPrev2(2), velPrev3(2), velPrev4(2), velExtrap(2);
    VectorXd  res(3), dp(2), Du(2), vel(2), velDot(2), force(2), res2(2), rStab(2), gradTvel(2), velPrev(2);
    MatrixXd  Dj(2, 3), grad(2,2), gradN(2,2), stress(2,2), gradPrev(2,2), matG(3,3), matJ(2,2), matJinv(2,2);

    elmDat = &(SolnData->ElemProp[elmType].data[0]);
    //matDat = &(SolidSolnData->MatlProp[matType].data[0]);

    rho = elmDat[4];
    mu  = elmDat[5];

    bool axsy = false;
    axsy = ((int)elmDat[2] == 1);

    af = SolnData->td(2);
    am = SolnData->td(3);
    acceFact = am*SolnData->td(9);
    muTaf = mu*af;
    dt = mpapTime.dt;

    volume = 0.5*(xNode[0]*(yNode[1]-yNode[2]) + xNode[1]*(yNode[2]-yNode[0]) + xNode[2]*(yNode[0]-yNode[1]));

    //cout << " volume = " << volume << endl;
    h2 = 4.0*volume/PI;
    h = sqrt(h2);

    stabParam = h2/(12.0*mu);
    tau[0] = elmDat[8]*stabParam;  // SUPG
    tau[1] = elmDat[9]*stabParam;  // PSPG
    tau[2] = elmDat[10]*stabParam; // LSIC

    nGP = elmDat[0];

    getGaussPointsTriangle(nGP, GeomData->gausspoints1, GeomData->gausspoints2, GeomData->gaussweights1);

    if(Klocal.rows() != nsize)
    {
      Klocal.resize(nsize, nsize);
      Flocal.resize(nsize);
    }

    Klocal.setZero();
    Flocal.setZero();
    
    //cout << " AAAAAAAAAA " << endl;
    //cout << nGP1 << '\t' << nGP2 << endl;

    for(gp1=0; gp1<nGP; gp1++)
    {
          param[0] = GeomData->gausspoints1[gp1];
          param[1] = GeomData->gausspoints2[gp1];

          GeomData->computeBasisFunctions2D(0, 1, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

          dvol = GeomData->gaussweights1[gp1] * Jac;
          
          LagrangeBasisFunsTria(degree, param[0], param[1], &N(0), dN_du[0], dN_du[1]);
          
          matJ.setZero();
          for(ii=0;ii<npElem;ii++)
          {
            matJ(0,0) +=  (xNode[ii] * dN_du[0][ii]) ;
            matJ(1,0) +=  (xNode[ii] * dN_du[1][ii]) ;
            matJ(0,1) +=  (yNode[ii] * dN_du[0][ii]) ;
            matJ(1,1) +=  (yNode[ii] * dN_du[1][ii]) ;
          }

          matJinv = matJ.inverse();

          matG.setZero();
          matG(0,0) = matJinv(0,0);
          matG(0,1) = matJinv(0,1);
          matG(1,0) = matJinv(1,0);
          matG(1,1) = matJinv(1,1);
          
          //printf("%12.6f \t %12.6f \t %12.6f \n", volume, h, Jac);
          //printf("%12.6f \t %12.6f \t %12.6f \t %12.6f \n", matB(0,0), matB(0,1), matB(1,0), matB(1,1));
          //printf("%12.6f \t %12.6f \t %12.6f \t %12.6f \n", matG(0,0), matG(0,1), matG(1,0), matG(1,1));

          matG = matG.transpose()*matG;

          vel(0) = computeValueCur(0, N);
          vel(1) = computeValueCur(1, N);

          grad(0,0) = computeValueCur(0, dN_dx);
          grad(0,1) = computeValueCur(0, dN_dy);
          grad(1,0) = computeValueCur(1, dN_dx);
          grad(1,1) = computeValueCur(1, dN_dy);

          velPrev(0) = computeValuePrev(0, N);
          velPrev(1) = computeValuePrev(1, N);

          velPrev2(0) = computeValuePrev2(0, N);
          velPrev2(1) = computeValuePrev2(1, N);

          velPrev3(0) = computeValuePrev3(0, N);
          velPrev3(1) = computeValuePrev3(1, N);

          velPrev4(0) = computeValuePrev4(0, N);
          velPrev4(1) = computeValuePrev4(1, N);

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
          //force(0) = analy.computeForce(0, xx, yy);
          //force(1) = analy.computeForce(1, xx, yy);
          //cout << force(0) << '\t' << force(1) << endl;

          //force(0) = computeForce(0, N);
          //force(1) = computeForce(1, N);

          velExtrap = velPrev;
          //velExtrap = 2.0*velPrev - velPrev2;
          //velExtrap = 3.0*velPrev - 3.0*velPrev2 + velPrev3;
          //velExtrap = 4.0*velPrev - 6.0*velPrev2 + 4.0*velPrev3 - velPrev4;

          gradTvel = grad*velExtrap ;

          res2(0) = rho*(velDot(0) + gradTvel(0) - force(0)) ;
          res2(1) = rho*(velDot(1) + gradTvel(1) - force(1)) ;

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

          velTemp(0) = velPrev(0);
          velTemp(1) = velPrev(1);
          velTemp(2) = 0.0;

          //evaluateStabParams_algo1(&velTemp(0), h, rho, mu, dt,  beta, tau);

          //evaluateStabParams_algo2(&velTemp(0), h, rho, mu, dt,  beta, tau);

          //evaluateStabParams_algo3(velTemp, matG, dt, rho, mu, CI, tau);

          //tau[0] *= elmDat[8];  // SUPG
          //tau[1] *= elmDat[9];  // PSPG
          //tau[2] *= elmDat[10]; // LSIC

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

            //Da = rho*(velPrev(0)*b1 + velPrev(1)*b2)*tau[0];
            Da = rho*(velExtrap(0)*b1 + velExtrap(1)*b2)*tau[0];

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

              // convection term - semi-implicit type A

              //Db = rho*(velPrev(0)*dN_dx(jj) + velPrev(1)*dN_dy(jj));
              Db = rho*(velExtrap(0)*dN_dx(jj) + velExtrap(1)*dN_dy(jj));

              Klocal(TI,   TJ)   += (b8*Db);
              Klocal(TIp1, TJp1) += (b8*Db);

              // pressure term
              Klocal(TI,   TJp2) -= (b1*N(jj));
              Klocal(TIp1, TJp2) -= (b2*N(jj));

              // continuity equation
              Klocal(TIp2, TJ)   += (b8*dN_dx(jj));
              Klocal(TIp2, TJp1) += (b8*dN_dy(jj));

              // SUPG and PSPG stabilisation terms

              Db *= af;

              Dj(0,0) = Db + fact2;
              Dj(0,1) = 0.0;
              Dj(0,2) = dN_dx(jj);
              Dj(1,0) = 0.0;
              Dj(1,1) = Db + fact2;
              Dj(1,2) = dN_dy(jj);

              if(axsy)
              {
                Dj(0,0) -= muTaf*(dN_dx(jj)/rad - N(jj)/rad/rad);
                Dj(1,1) -= muTaf*(dN_dx(jj)/rad);
              }

              // SUPG
              Klocal(TI, TJ)     += Da*Dj(0,0);
              Klocal(TI, TJp1)   += Da*Dj(0,1);
              Klocal(TI, TJp2)   += Da*Dj(0,2);

              Klocal(TIp1, TJ)   += Da*Dj(1,0);
              Klocal(TIp1, TJp1) += Da*Dj(1,1);
              Klocal(TIp1, TJp2) += Da*Dj(1,2);

              //Klocal(TI,   TJ)   += ( (tau[0]*af) * b1 * rStab(0) * N(jj) );
              //Klocal(TI,   TJp1) += ( (tau[0]*af) * b2 * rStab(0) * N(jj) );
              //Klocal(TIp1, TJ)   += ( (tau[0]*af) * b1 * rStab(1) * N(jj) );
              //Klocal(TIp1, TJp1) += ( (tau[0]*af) * b2 * rStab(1) * N(jj) );

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
*/



/*
int LagrangeElem2DNavierStokesTria3Node::calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal)
{
//   char fct[] = "LagrangeElem2DNavierStokes4Node::calcStiffnessAndResidual";
//   computerTime.go(fct);

  // stabilised formulation
  // semi-implicit scheme - type B

    int ii, jj, gp1, gp2, TI, TIp1, TIp2, count, TJ, TJp1, TJp2;

    double  uu, vv, Jac, dvol, b1, b2, b3, b4, b5, b6, b7, b8, xx, yy, acceFact, rho, mu, CI=1.0;
    double  pres, Da, Db, af, am, d1, c1, muTaf, rad, urdr, urdr2, tau[3], volume;
    double  fact, fact1, fact2, bb1, bb2, param[2], dt, h, h2, stabParam, *elmDat, x[3], y[3], beta[6];
    double  dN_du[2][3];

    beta[0] = 1.0;     beta[1] = 1.0/3.0;
    beta[2] = 30.0;    beta[3] = 0.1;
    beta[4] = 1.0;     beta[5] = 1.0;

    for(ii=0;ii<npElem;ii++)
    {
      x[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      y[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    VectorXd  N(nlbf), dN_dx(nlbf), dN_dy(nlbf), velTemp(3);
    VectorXd  res(3), dp(2), Du(2), vel(2), velDot(2), force(2), res2(2), rStab(2), gradTvel(2), velPrev(2);
    MatrixXd  Dj(2, 3), grad(2,2), gradN(2,2), stress(2,2), gradPrev(2,2), matG(3,3), matJ(2,2), matJinv(2,2);

    elmDat = &(SolnData->ElemProp[elmType].data[0]);
    //matDat = &(SolidSolnData->MatlProp.data[0]);

    rho = elmDat[4];
    mu  = elmDat[5];

    bool axsy = false;
    axsy = ((int)elmDat[2] == 1);

    af = SolnData->td(2);
    am = SolnData->td(3);
    acceFact = am*SolnData->td(9);
    muTaf = mu*af;
    dt = mpapTime.dt;

    volume = 0.5*(x[0]*(y[1]-y[2]) + x[1]*(y[2]-y[0]) + x[2]*(y[0]-y[1]));

    //cout << " volume = " << volume << endl;
    h2 = 4.0*volume/PI;
    h = sqrt(h2);

    stabParam = h2/(12.0*mu);
    tau[0] = elmDat[8]*stabParam;  // SUPG
    tau[1] = elmDat[9]*stabParam;  // PSPG
    tau[2] = elmDat[10]*stabParam; // LSIC

    nGP = elmDat[0];

    getGaussPointsTriangle(nGP, GeomData->gausspoints1, GeomData->gausspoints2, GeomData->gaussweights1);

    if(Klocal.rows() != nsize)
    {
      Klocal.resize(nsize, nsize);
      Flocal.resize(nsize);
    }
    Klocal.setZero();
    Flocal.setZero();

    //cout << " AAAAAAAAAA " << endl;
    //cout << nGP << '\t' << nGP << endl;

    for(gp1=0;gp1<nGP;gp1++)
    {
          param[0] = GeomData->gausspoints1[gp1];
          param[1] = GeomData->gausspoints2[gp1];

          GeomData->computeBasisFunctions2D(0, 1, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

          dvol = GeomData->gaussweights1[gp1] * Jac;
          
          LagrangeBasisFunsTria(degree, param[0], param[1], &N(0), dN_du[0], dN_du[1]);
          
          //matJ.setZero();
          //for(ii=0;ii<npElem;ii++)
          //{
            //matJ(0,0) +=  (x[ii] * dN_du[0][ii]) ;
            //matJ(1,0) +=  (x[ii] * dN_du[1][ii]) ;
            //matJ(0,1) +=  (y[ii] * dN_du[0][ii]) ;
            //matJ(1,1) +=  (y[ii] * dN_du[1][ii]) ;
          //}

          //matJinv = matJ.inverse();

          //matG.setZero();
          //matG(0,0) = matJinv(0,0);
          //matG(0,1) = matJinv(0,1);
          //matG(1,0) = matJinv(1,0);
          //matG(1,1) = matJinv(1,1);
          
          //printf("%12.6f \t %12.6f \t %12.6f \n", volume, h, Jac);
          //printf("%12.6f \t %12.6f \t %12.6f \t %12.6f \n", matB(0,0), matB(0,1), matB(1,0), matB(1,1));
          //printf("%12.6f \t %12.6f \t %12.6f \t %12.6f \n", matG(0,0), matG(0,1), matG(1,0), matG(1,1));

          //matG = matG.transpose()*matG;

          vel(0) = computeValueCur(0, N);
          vel(1) = computeValueCur(1, N);

          grad(0,0) = computeValueCur(0, dN_dx);
          grad(0,1) = computeValueCur(0, dN_dy);
          grad(1,0) = computeValueCur(1, dN_dx);
          grad(1,1) = computeValueCur(1, dN_dy);

          velPrev(0) = computeValuePrev(0, N);
          velPrev(1) = computeValuePrev(1, N);

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
          //force(0) = analy.computeForce(0, xx, yy);
          //force(1) = analy.computeForce(1, xx, yy);
          //cout << force(0) << '\t' << force(1) << endl;

          //force(0) = computeForce(0, N);
          //force(1) = computeForce(1, N);


          gradTvel = gradPrev*vel + grad*velPrev - gradPrev*velPrev;

          res2(0) = rho*(velDot(0) + gradTvel(0) - force(0)) ;
          res2(1) = rho*(velDot(1) + gradTvel(1) - force(1)) ;

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

          velTemp(0) = velPrev(0);
          velTemp(1) = velPrev(1);
          velTemp(2) = 0.0;

          //evaluateStabParams_algo1(&velTemp(0), h, rho, mu, dt,  beta, tau);

          //evaluateStabParams_algo2(&velTemp(0), h, rho, mu, dt,  beta, tau);

          //evaluateStabParams_algo3(velTemp, matG, dt, rho, mu, CI, tau);

          //tau[0] *= elmDat[8];  // SUPG
          //tau[1] *= elmDat[9];  // PSPG
          //tau[2] *= elmDat[10]; // LSIC

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
            //Da = rho*(vel(0)*b1 + vel(1)*b2)*tau[0];

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

              // convection term

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

              gradN *= af;

              Dj(0,0) = gradN(0,0) + fact2;
              Dj(0,1) = gradN(0,1);
              Dj(0,2) = dN_dx(jj);
              Dj(1,0) = gradN(1,0);
              Dj(1,1) = gradN(1,1) + fact2;
              Dj(1,2) = dN_dy(jj);

              if(axsy)
              {
                Dj(0,0) -= muTaf*(dN_dx(jj)/rad - N(jj)/rad/rad);
                Dj(1,1) -= muTaf*(dN_dx(jj)/rad);
              }

              // SUPG
              Klocal(TI, TJ)     += Da*Dj(0,0);
              Klocal(TI, TJp1)   += Da*Dj(0,1);
              Klocal(TI, TJp2)   += Da*Dj(0,2);

              Klocal(TIp1, TJ)   += Da*Dj(1,0);
              Klocal(TIp1, TJp1) += Da*Dj(1,1);
              Klocal(TIp1, TJp2) += Da*Dj(1,2);

              //Klocal(TI,   TJ)   += ( (tau[0]*af) * b1 * rStab(0) * N(jj) );
              //Klocal(TI,   TJp1) += ( (tau[0]*af) * b2 * rStab(0) * N(jj) );
              //Klocal(TIp1, TJ)   += ( (tau[0]*af) * b1 * rStab(1) * N(jj) );
              //Klocal(TIp1, TJp1) += ( (tau[0]*af) * b2 * rStab(1) * N(jj) );

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
*/




int LagrangeElem2DNavierStokesTria3Node::calcInternalForces()
{
  return 0;
}



void LagrangeElem2DNavierStokesTria3Node::discreteContourplot(int vartype, int varindex, int index, int nCol, double umin, double umax)
{
  return;
}


void LagrangeElem2DNavierStokesTria3Node::projectToKnots(bool extrapolateFlag, int vartype, int varindex, int index)
{
  return;
}


void LagrangeElem2DNavierStokesTria3Node::projectStress(int varindex, double* outval)
{
  return;
}



void LagrangeElem2DNavierStokesTria3Node::projectStrain(int vartype, int varindex, double* outval)
{
  return;
}



void LagrangeElem2DNavierStokesTria3Node::projectIntVar(int index, double* outval)
{
  return;
}


int LagrangeElem2DNavierStokesTria3Node::calcOutput(double u1, double v1)
{
  return 0;
}



void LagrangeElem2DNavierStokesTria3Node::toPostprocess(int vartype, int varindex, int type, SparseMatrixXd&  coeffMat, VectorXd& rhsVec)
{
  return;
}

