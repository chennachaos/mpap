
#include <math.h>
#include "Debug.h"
//#include "Plot.h"
//#include "FunctionsElement.h"
#include "MpapTime.h"
#include "NurbsElem2DStructMixed2fieldStabilised.h"
#include "NurbsShapeFunctions.h"
#include <assert.h>
#include "ComputerTime.h"

#include "util.h"

#include "ExactSolutionsElasticity.h"

using namespace std;

//extern Plot     plot;
extern MpapTime mpapTime;
extern ComputerTime       computerTime;



NurbsElem2DStructMixed2fieldStabilised::NurbsElem2DStructMixed2fieldStabilised()
{
  if (debug) cout << " constructor NurbsElem2DStructMixed2fieldStabilised\n\n";

//   cout << " constructor NurbsElem2DStructMixed2fieldStabilised\n\n";

  calcExtraMatrices = true;
}



NurbsElem2DStructMixed2fieldStabilised::~NurbsElem2DStructMixed2fieldStabilised()
{
  if (debug) cout << " destructor NurbsElem2DStructMixed2fieldStabilised\n\n";

  // cout << " destructor NurbsElem2DStructMixed2fieldStabilised\n\n";

}

int NurbsElem2DStructMixed2fieldStabilised::calcStiffnessMatrix(double dt)
{

  return 0;
}



void NurbsElem2DStructMixed2fieldStabilised::contourplot(int index, int nCol, double umin, double umax)
{
  return;
}

int NurbsElem2DStructMixed2fieldStabilised::calcMassMatrix(int lumpInd, double dt)
{


  return 0;
}


int NurbsElem2DStructMixed2fieldStabilised::calcStiffnessAndResidual()
{
  if(axsy)
  {
    //cout << " AXISYMMETRIC case is not supported by NurbsElem2DStructMixed2fieldStabilised  " << endl;

    if(finite)
      NurbsElem2DStructMixed2fieldStabilised::calcStiffnessAndResidualAxsyFS();
    else
      NurbsElem2DStructMixed2fieldStabilised::calcStiffnessAndResidualAxsySS();
  }
  else
  {
    if(finite)
      NurbsElem2DStructMixed2fieldStabilised::calcStiffnessAndResidual2();
    else
      NurbsElem2DStructMixed2fieldStabilised::calcStiffnessAndResidual1();
  }

  return 0;
}


/*
int NurbsElem2DStructMixed2fieldStabilised::calcStiffnessAndResidual1()
{
  // stabilized formulation for small strains

  int  err, isw, count, count1, ll, ii, jj, TI, TIp1, TIp2, TJ, TJp1, TJp2, index, gp1, gp2, mm;

  double  F[4], detF, F33, fact, dvol0, dt, Jac, dummy, pres, volstr, force[2];
  double  cc[4][4], stre[4], bc[3][3], Idev[4][4], cctmp[4][4];
  double  N[nlbf], dN_dx[nlbf], dN_dy[nlbf], d2N_dx2[nlbf], d2N_dy2[nlbf], d2N_dxy[nlbf], d2N_dyx[nlbf];
  double  eps, BULK = matDat[0], tau, mu, dp[2], h2, cc1, cc2, cc3, alpha, r2D3=2.0/3.0;
  double  r2mu, bb1, bb2, bb3, bb4, bb5, ci;
  
  MatrixXd  D2u(2,3), Ka(3,3), Kb(3,3), Kstab(3,3);
  VectorXd  rStab(3);
  MatrixXd  D(nsize,3);

  mu = matDat[1];
  
  r2mu = 2.0*mu;
  
  h2  = volume()*4.0/PI;
//  h2  = sqrt(volume());

  alpha = 1.0;
//  alpha = 1.0/3.0;
//  alpha = 1.0/10.0;
//  alpha = 1.0/12.0;
//  alpha = 1.0/20.0;
//  alpha = 1.0/40.0;
//  alpha = 1.0/50.0;
//  alpha = 1.0/100.0;
//  alpha *= 1.0/surf0->p/surf0->q;

  tau = elmDat[6]*alpha*h2/mu;
  
  //cout << " tau = " << volume() << '\t' << tau << endl;

  Idev2D(Idev);

  double *gaussweights = &(surf0->gaussweights[0]);

  double  *values1 = &(surf1->Values[0][0]);
  double  *values2 = &(surf1->Values[1][0]);
  double  *values3 = &(surf1->Values[2][0]);

  int *tt = &(surf0->IEN[elenum][0]);

  eps = 1.0/BULK;

  Klocal.setZero();
  Flocal.setZero();

  count = 1;   ll = 0;   err = 0;   isw = 3;
  dt = mpapTime.dt;

  count1 = 0;
  for(gp2=0;gp2<nGP2;gp2++)
  {
      for(gp1=0;gp1<nGP1;gp1++)
      {
          index = count1*2;

          surf0->ShapeFunDerivatives2(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, d2N_dx2, d2N_dy2, d2N_dxy, d2N_dyx, Jac);
          //cout << " bbbbbbbbbbbbb " << endl;
          dvol0 = Jac * gaussweights[count1] * JacMultFact;

          //cout << " bbbbbbbbbbbbb " << endl;
          //surf1->deformationGradient(startindex[0], startindex[1], 1, dN_dx, dN_dy, F, detF);

          F[0] = F[1] = F[2] = F[3] = 0.0;
          pres = 0.0;
          dp[0] = dp[1] = 0.0;
          D2u.setZero();
          D.setZero();
          for(ii=0;ii<nlbf;ii++)
          {
            index = tt[ii];

            bb1 = values1[index];
            bb2 = values2[index];
            bb3 = values3[index];

            pres  += N[ii]*bb3;

            dp[0] += dN_dx[ii]*bb3;
            dp[1] += dN_dy[ii]*bb3;

            F[0] += bb1*dN_dx[ii];
            F[2] += bb1*dN_dy[ii];
            F[1] += bb2*dN_dx[ii];
            F[3] += bb2*dN_dy[ii];

            D2u(0,0) += bb1*d2N_dx2[ii];
            D2u(0,1) += bb1*d2N_dxy[ii];
            D2u(0,2) += bb1*d2N_dy2[ii];

            D2u(1,0) += bb2*d2N_dx2[ii];
            D2u(1,1) += bb2*d2N_dxy[ii];
            D2u(1,2) += bb2*d2N_dy2[ii];
            
            TI   =  3*ii;
            TIp1 =  TI+1;
            TIp2 =  TI+2;

            bb1 = dN_dx[ii];
            bb2 = dN_dy[ii];

            fact = d2N_dx2[ii]+d2N_dy2[ii];

            ci = mu * fact;

            //D(TI,0)   = ci;
            //D(TIp1,0) = 0.0;
            //D(TIp2,0) = bb1;

            //D(TI,1)   = 0.0;
            //D(TIp1,1) = ci;
            //D(TIp2,1) = bb2;

            //D(TI,2)   =  bb1;
            //D(TIp1,2) =  bb2;
            //D(TIp2,2) = -N[ii]*eps;

            D(TI,   0) = r2mu*(r2D3*d2N_dx2[ii]+d2N_dy2[ii]);
            D(TIp1, 0) = r2mu*(r2D3*d2N_dxy[ii]);
            D(TIp2, 0) = dN_dx[ii];

            D(TI,   1) = r2mu*(r2D3*d2N_dxy[ii]);
            D(TIp1, 1) = r2mu*(d2N_dx2[ii]+r2D3*d2N_dy2[ii]);
            D(TIp2, 1) = dN_dy[ii];

            //D(TI,   2) = dN_dx[ii];
            //D(TIp1, 2) = dN_dy[ii];
            //D(TIp2, 2) = -eps*N[ii];

            D(TI,   2) = 0.0;
            D(TIp1, 2) = 0.0;
            D(TIp2, 2) = 0.0;

            //printf(" d2N ");        printf("\t%12.8f\t%12.8f\t%12.8f\n", d2N_dx2[ii], d2N_dxy[ii], d2N_dy2[ii]);
          }
          F[0] += 1.0;
          F[3] += 1.0;

          //cout << " cccccccccccccc " << endl;

          volstr = (F[0]+F[3]-2.0);

          F33 = 1.0;

          //cout << " ddddddddddddddd " << endl;

          matlib2d_(matDat, F, &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);
          if(err !=0)           return 1;

          dummy = pres - (stre[0]+stre[1]+stre[2])/3.0 ;

          stre[0] += dummy;
          stre[1] += dummy;
          stre[2] += dummy;

          //printf(" stresses ");        printf("\t%12.8f\t%12.8f\t%12.8f\t%12.8f\n\n", stre[0], stre[1], stre[2], pres);

          for(ii=0;ii<4;ii++)
          {
             for(jj=0;jj<4;jj++)
             {
                cctmp[ii][jj] = 0.0;
                for(mm=0;mm<4;mm++)
                   cctmp[ii][jj] += Idev[ii][mm]*cc[mm][jj];
             }
          }

          for(ii=0;ii<4;ii++)
          {
             for(jj=0;jj<4;jj++)
             {
               cc[ii][jj] = 0.0;
               for(mm=0;mm<4;mm++)
                  cc[ii][jj] += cctmp[ii][mm]*Idev[mm][jj];
             }
          }

          //==============================================
          // CALCULATE TANGENT STIFFNESS and RESIDUAL
          //==============================================

          //force[0] = bforce[0] * rho0 ;
          //force[1] = bforce[1] * rho0 ;
          force[0] = 0.0;
          force[1] = 0.0;

          fact = volstr - pres*eps;

          for(ii=0;ii<nlbf;ii++)
          {
            bb1 = dN_dx[ii] * dvol0;
            bb2 = dN_dy[ii] * dvol0;
            bb3 = N[ii] * dvol0;

            bc[0][0] = (bb1 * cc[0][0] + bb2 * cc[3][0]);
            bc[0][1] = (bb1 * cc[0][1] + bb2 * cc[3][1]);
            bc[0][2] = (bb1 * cc[0][3] + bb2 * cc[3][3]);

            bc[1][0] = (bb2 * cc[1][0] + bb1 * cc[3][0]);
            bc[1][1] = (bb2 * cc[1][1] + bb1 * cc[3][1]);
            bc[1][2] = (bb2 * cc[1][3] + bb1 * cc[3][3]);

            TI   = 3*ii;
            TIp1 = TI+1;
            TIp2 = TI+2;

            Flocal[TI]   += (bb3*force[0] - bb1*stre[0] - bb2*stre[3]) ;
            Flocal[TIp1] += (bb3*force[1] - bb1*stre[3] - bb2*stre[1]) ;
            Flocal[TIp2] -= (bb3*fact);

            for(jj=0;jj<nlbf;jj++)
            {
                TJ   = 3*jj;
                TJp1 = TJ+1;
                TJp2 = TJ+2;

                cc1 = dN_dx[jj];
                cc2 = dN_dy[jj];
                cc3 = N[jj];

                Klocal(TI,   TJ)    +=  (bc[0][0] * cc1 + bc[0][2] * cc2) ;
                Klocal(TI,   TJp1)  +=  (bc[0][1] * cc2 + bc[0][2] * cc1) ;
                Klocal(TIp1, TJ)    +=  (bc[1][0] * cc1 + bc[1][2] * cc2) ;
                Klocal(TIp1, TJp1)  +=  (bc[1][1] * cc2 + bc[1][2] * cc1) ;

                Klocal(TI,   TJp2)  +=  bb1 * cc3 ;
                Klocal(TIp1, TJp2)  +=  bb2 * cc3 ;

                Klocal(TIp2, TJ)    +=  bb3 * cc1 ;
                Klocal(TIp2, TJp1)  +=  bb3 * cc2 ;

                Klocal(TIp2, TJp2)  -= eps*(bb3*cc3);
            }
          }


          //printf(" %14.8f \t%14.8f \t%14.8f \t%14.8f \t%14.8f \t%14.8f \t%14.8f \t%14.8f\n\n", ff, Du(0), Du(1), dp(0), dp(1), f(0), f(1), pres);

          rStab(0) = force[0] - (r2mu*(r2D3*D2u(0,0) + D2u(0,2) + r2D3*D2u(1,1)) + dp[0]);
          rStab(1) = force[1] - (r2mu*(r2D3*D2u(0,1) + D2u(1,0) + r2D3*D2u(1,2)) + dp[1]);
          //rStab(2) = 0.0 - (volstr - pres*eps);
          rStab(2) = 0.0 ;

          //rStab(0) = force[0] - (mu*(D2u(0,0) + D2u(0,2)) + dp[0]);
          //rStab(1) = force[1] - (mu*(D2u(1,0) + D2u(1,2)) + dp[1]);
          //rStab(2) = 0.0 - (volstr - pres*eps);
          //rStab(2) = 0.0;

          dvol0 *= tau;

          for(ii=0;ii<nsize;ii++)
          {
            Flocal[ii] -= dvol0*(D(ii,0)*rStab(0)+D(ii,1)*rStab(1)+D(ii,2)*rStab(2));

            for(jj=0;jj<nsize;jj++)
               Klocal(ii, jj)  -=  dvol0*(D(ii,0)*D(jj,0)+D(ii,1)*D(jj,1) + D(ii,2)*D(jj,2)) ;
          }

          count++;
          count1++;
          ll += nivGP;
     }
   }

//  printStiffnessMatrix();
//  printf("\n\n");
//  printForceVector();

  return 0;
}
*/



int NurbsElem2DStructMixed2fieldStabilised::calcStiffnessAndResidual1()
{
  // PSPG stabilized formulation for small strains

  int  err, isw, count, count1, ll, ii, jj, kk, TI, TIp1, TIp2, TJ, TJp1, TJp2, index, gp1, gp2, mm;

  double  F[4], detF, F33, fact, fact2, dvol0, dt, Jac, dummy, pres, volstr, force[2];
  double  cc[4][4], stre[4], bc[3][3], Idev[4][4], cctmp[4][4];
  double  N[nlbf], dN_dx[nlbf], dN_dy[nlbf], d2N_dx2[nlbf], d2N_dy2[nlbf], d2N_dxy[nlbf], d2N_dyx[nlbf];
  double  eps, BULK = matDat[0], tau, mu, dp[2], h2, cc1, cc2, cc3, alpha, r2D3=2.0/3.0;
  double  r2mu, bb1, bb2, bb3, bb4, bb5, ci;
  double  veloCur[2], acceCur[2], acceFact;
  
  MatrixXd  D2u(2,3), Ka(3,3), Kb(3,3), Kstab(3,3), Dj(2,3);
  VectorXd  rStab(2);
  MatrixXd  D(nsize,3);

  BULK = matDat[0];
  mu   = matDat[1];
  eps  = 1.0/BULK;

  r2mu = 2.0*mu;
  
  h2  = volume()*4.0/PI;
//  h2  = sqrt(volume());

//  alpha = 1.0;
//  alpha = 1.0/3.0;
//  alpha = 1.0/10.0;
  alpha = 1.0/12.0;
//  alpha = 1.0/20.0;
//  alpha = 1.0/40.0;
//  alpha = 1.0/50.0;
//  alpha = 1.0/100.0;
  alpha *= 1.0/surf0->p/surf0->q;

  tau = elmDat[8]*alpha*h2/mu;

  //cout << " tau = " << volume() << '\t' << tau << endl;
  //cout << " nGP = " << nGP1 << '\t' << nGP2 << endl;

  Idev2D(Idev);

  int *tt = &(surf0->IEN[elenum][0]);

  double *gaussweights1 = &(surf0->gaussweights1[0]);
  double *gaussweights2 = &(surf0->gaussweights2[0]);

  double  *valuesCur1 = &(surf1->ValuesCur[0][0]);
  double  *valuesCur2 = &(surf1->ValuesCur[1][0]);
  double  *valuesCur3 = &(surf1->ValuesCur[2][0]);

  double  *valuesDotCur1 = &(surf1->ValuesDotCur[0][0]);
  double  *valuesDotCur2 = &(surf1->ValuesDotCur[1][0]);

  double  *valuesDotDotCur1 = &(surf1->ValuesDotDotCur[0][0]);
  double  *valuesDotDotCur2 = &(surf1->ValuesDotDotCur[1][0]);


  double rho0 = elmDat[7];
  double  af = surf1->td(2);
  double  d1 = surf1->td(5);
  double  aa = surf1->td(10);


  Klocal.setZero();
  Flocal.setZero();

  count = 1;   ll = 0;   err = 0;   isw = 3;
  dt = mpapTime.dt;

  count1 = 0;
  for(gp2=0;gp2<nGP2;gp2++)
  {
      for(gp1=0;gp1<nGP1;gp1++)
      {
          index = count1*2;

          //surf0->ShapeFunDerivatives2(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, d2N_dx2, d2N_dy2, d2N_dxy, d2N_dyx, Jac);
          surf0->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), &N[0], &dN_dx[0], &dN_dy[0], Jac);
          //cout << " bbbbbbbbbbbbb " << endl;
          fact = gaussweights2[gp2] * gaussweights1[gp1] * thick * JacMultFact;

          dvol0 = Jac * fact;

          //cout << " bbbbbbbbbbbbb " << endl;
          //surf1->deformationGradient(startindex[0], startindex[1], 1, dN_dx, dN_dy, F, detF);

          veloCur[0] = veloCur[1] = 0.0;
          acceCur[0] = acceCur[1] = 0.0;

          F[0] = F[1] = F[2] = F[3] = 0.0;
          pres = 0.0;
          dp[0] = dp[1] = 0.0;
          D2u.setZero();
          D.setZero();
          for(ii=0;ii<nlbf;ii++)
          {
            kk = tt[ii];

            bb1 = valuesCur1[kk];
            bb2 = valuesCur2[kk];
            bb3 = valuesCur3[kk];

            pres  += N[ii]*bb3;

            dp[0] += dN_dx[ii]*bb3;
            dp[1] += dN_dy[ii]*bb3;

            F[0] += bb1*dN_dx[ii];
            F[2] += bb1*dN_dy[ii];
            F[1] += bb2*dN_dx[ii];
            F[3] += bb2*dN_dy[ii];

            D2u(0,0) += bb1*d2N_dx2[ii];
            D2u(0,1) += bb1*d2N_dxy[ii];
            D2u(0,2) += bb1*d2N_dy2[ii];

            D2u(1,0) += bb2*d2N_dx2[ii];
            D2u(1,1) += bb2*d2N_dxy[ii];
            D2u(1,2) += bb2*d2N_dy2[ii];
            
            TI   =  3*ii;
            TIp1 =  TI+1;
            TIp2 =  TI+2;

            bb1 = dN_dx[ii];
            bb2 = dN_dy[ii];

            fact = d2N_dx2[ii]+d2N_dy2[ii];

            ci = mu * fact;

            //D(TI,2)   =  bb1;
            //D(TIp1,2) =  bb2;
            //D(TIp2,2) = -N[ii]*eps;

            D(TI,   0) = r2mu*(r2D3*d2N_dx2[ii]+d2N_dy2[ii]);
            D(TIp1, 0) = r2mu*(r2D3*d2N_dxy[ii]);
            D(TIp2, 0) = dN_dx[ii];

            D(TI,   1) = r2mu*(r2D3*d2N_dxy[ii]);
            D(TIp1, 1) = r2mu*(d2N_dx2[ii]+r2D3*d2N_dy2[ii]);
            D(TIp2, 1) = dN_dy[ii];

            veloCur[0] += valuesDotCur1[kk]*N[ii];
            veloCur[1] += valuesDotCur2[kk]*N[ii];

            acceCur[0] += valuesDotDotCur1[kk]*N[ii];
            acceCur[1] += valuesDotDotCur2[kk]*N[ii];

            //printf(" d2N ");        printf("\t%12.8f\t%12.8f\t%12.8f\n", d2N_dx2[ii], d2N_dxy[ii], d2N_dy2[ii]);
          }
          F[0] += 1.0;
          F[3] += 1.0;

          //cout << " cccccccccccccc " << endl;

          volstr = (F[0]+F[3]-2.0);

          F33 = 1.0;

          //cout << " ddddddddddddddd " << endl;

          matlib2d_(matDat, F, &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);
          if(err !=0)           return 1;

          dummy = pres - (stre[0]+stre[1]+stre[2])/3.0 ;

          stre[0] += dummy;
          stre[1] += dummy;
          stre[2] += dummy;

          //printf(" stresses ... \t%12.8f\t%12.8f\t%12.8f\t%12.8f\n\n", stre[0], stre[1], stre[2], pres);

          for(ii=0;ii<4;ii++)
          {
             for(jj=0;jj<4;jj++)
             {
                cctmp[ii][jj] = 0.0;
                for(mm=0;mm<4;mm++)
                   cctmp[ii][jj] += Idev[ii][mm]*cc[mm][jj];
             }
          }

          for(ii=0;ii<4;ii++)
          {
             for(jj=0;jj<4;jj++)
             {
               cc[ii][jj] = 0.0;
               for(mm=0;mm<4;mm++)
                  cc[ii][jj] += cctmp[ii][mm]*Idev[mm][jj];
             }
          }

          //==============================================
          // CALCULATE TANGENT STIFFNESS and RESIDUAL
          //==============================================

          //force[0] = bforce[0] * rho0 ;
          //force[1] = bforce[1] * rho0 ;
          force[0] = 0.0;
          force[1] = 0.0;
          force[0] += rho0*acceCur[0];
          force[1] += rho0*acceCur[1];

          rStab[0] = rho0*acceCur[0] - dp[0];
          rStab[1] = rho0*acceCur[1] - dp[1];

          volstr = volstr - pres*eps;

          for(ii=0;ii<nlbf;ii++)
          {
            bb1 = dN_dx[ii] * dvol0;
            bb2 = dN_dy[ii] * dvol0;
            bb3 = N[ii] * dvol0;

            bc[0][0] = (bb1 * cc[0][0] + bb2 * cc[3][0]);
            bc[0][1] = (bb1 * cc[0][1] + bb2 * cc[3][1]);
            bc[0][2] = (bb1 * cc[0][3] + bb2 * cc[3][3]);

            bc[1][0] = (bb2 * cc[1][0] + bb1 * cc[3][0]);
            bc[1][1] = (bb2 * cc[1][1] + bb1 * cc[3][1]);
            bc[1][2] = (bb2 * cc[1][3] + bb1 * cc[3][3]);

            TI   = 3*ii;
            TIp1 = TI+1;
            TIp2 = TI+2;

            Flocal[TI]   -= (bb3*force[0] + bb1*stre[0] + bb2*stre[3]) ;
            Flocal[TIp1] -= (bb3*force[1] + bb1*stre[3] + bb2*stre[1]) ;
            Flocal[TIp2] -= (bb3*volstr);

            // PSPG stabilization
            Flocal[TIp2] -= tau*(bb1*rStab[0] + bb2*rStab[1]);

            acceFact = d1*rho0*bb3;

            for(jj=0;jj<nlbf;jj++)
            {
                TJ   = 3*jj;
                TJp1 = TJ+1;
                TJp2 = TJ+2;

                cc1 = dN_dx[jj];
                cc2 = dN_dy[jj];
                cc3 = N[jj];

                fact2 = acceFact*N[jj];

                Klocal(TI,   TJ)    += fact2 ;
                Klocal(TIp1, TJp1)  += fact2 ;

                Klocal(TI,   TJ)    +=  af*(bc[0][0] * cc1 + bc[0][2] * cc2) ;
                Klocal(TI,   TJp1)  +=  af*(bc[0][1] * cc2 + bc[0][2] * cc1) ;
                Klocal(TI,   TJp2)  +=  af*(bb1 * cc3) ;

                Klocal(TIp1, TJ)    +=  af*(bc[1][0] * cc1 + bc[1][2] * cc2) ;
                Klocal(TIp1, TJp1)  +=  af*(bc[1][1] * cc2 + bc[1][2] * cc1) ;
                Klocal(TIp1, TJp2)  +=  af*(bb2 * cc3) ;

                Klocal(TIp2, TJ)    +=  af*(bb3 * cc1) ;
                Klocal(TIp2, TJp1)  +=  af*(bb3 * cc2) ;
                Klocal(TIp2, TJp2)  -=  af*(bb3 * cc3)*eps;

                // PSPG stabilization

                //fact2 -= ( muTaf*d2N(jj) );
                //fact2 = 0.0;
                fact2 = d1*rho0*cc3;

                Dj(0,0) = fact2;
                Dj(0,1) = 0.0;
                Dj(0,2) = -af*dN_dx[jj];
                Dj(1,0) = 0.0;
                Dj(1,1) = fact2;
                Dj(1,2) = -af*dN_dy[jj];

                Klocal(TIp2, TJ)   += (bb1*Dj(0,0) + bb2*Dj(1,0))*tau;
                Klocal(TIp2, TJp1) += (bb1*Dj(0,1) + bb2*Dj(1,1))*tau;
                Klocal(TIp2, TJp2) += (bb1*Dj(0,2) + bb2*Dj(1,2))*tau;
            }
          }

          //printf(" %14.8f \t%14.8f \t%14.8f \t%14.8f \t%14.8f \t%14.8f \t%14.8f \t%14.8f\n\n", ff, Du(0), Du(1), dp(0), dp(1), f(0), f(1), pres);

          //rStab(0) = force[0] - (r2mu*(r2D3*D2u(0,0) + D2u(0,2) + r2D3*D2u(1,1)) + dp[0]);
          //rStab(1) = force[1] - (r2mu*(r2D3*D2u(0,1) + D2u(1,0) + r2D3*D2u(1,2)) + dp[1]);
          //rStab(2) = 0.0 - (volstr - pres*eps);
          //rStab(2) = 0.0 ;

          //rStab(0) = force[0] - (mu*(D2u(0,0) + D2u(0,2)) + dp[0]);
          //rStab(1) = force[1] - (mu*(D2u(1,0) + D2u(1,2)) + dp[1]);
          //rStab(2) = 0.0 - (volstr - pres*eps);
          //rStab(2) = 0.0;

          count++;
          count1++;
          ll += nivGP;
     }
   }

  //printStiffnessMatrix();   printf("\n\n");   printForceVector();

  return 0;
}




int NurbsElem2DStructMixed2fieldStabilised::toComputeInfSupCondition()
{
  return 0;
}





int NurbsElem2DStructMixed2fieldStabilised::calcStiffnessAndResidual2()
{
  // PSPG stabilized formulation
  
  int  err, isw, count, count1, ll, ii, jj, kk, TI, TIp1, TIp2, TJ, TJp1, TJp2, index, gp1, gp2, mm;

  double  F[4], detF, F33, dvol0, dvol, dt, Jac, dummy, pres, pbar, r2d3 = 2.0/3.0;
  double  stre[4], strdev[4], bc[2][3], dN_dx[nlbf], dN_dy[nlbf], N[nlbf];
  double  fact, fact1, fact2, fact3, fact4;
  double  Idev[4][4], D11[4][4], cctmp[4][4], force[2], d2N_dx2[nlbf], d2N_dy2[nlbf], d2N_dxy[nlbf], d2N_dyx[nlbf];
  double  tau, dp[2], h, h2, alpha, sig[2];
  double  bb1, bb2, bb3, bb4, cc1, cc2, cc3, volstr, veloCur[2], acceCur[2], rStab[2], acceFact;

  MatrixXd  D2u(2,3), Dj(2,3);

  double rho0 = elmDat[7];
  double BULK = matDat[0];
  double mu   = matDat[1];
  double eps  = 1.0/BULK;

  h2  = volume()*4.0/PI;
  h = sqrt(h2);
  
  alpha = 1.0;
  alpha = 1.0/10.0;
  alpha = 1.0/20.0;
  alpha = 1.0/40.0;
  alpha = 1.0/12.0;
//  alpha *= 1.0/surf0->p/surf0->q;

  tau = elmDat[8]*alpha*h2/mu;

//  tau = elmDat[8]*h2/(2.0*sqrt(BULK+4.0*mu/3.0));

  //cout << " tau = " << volume() << '\t' << tau << endl;

  int *tt = &(surf0->IEN[elenum][0]);

  double *gaussweights1 = &(surf0->gaussweights1[0]);
  double *gaussweights2 = &(surf0->gaussweights2[0]);

  double  *valuesCur1 = &(surf1->ValuesCur[0][0]);
  double  *valuesCur2 = &(surf1->ValuesCur[1][0]);
  double  *valuesCur3 = &(surf1->ValuesCur[2][0]);

  double  *valuesDotCur1 = &(surf1->ValuesDotCur[0][0]);
  double  *valuesDotCur2 = &(surf1->ValuesDotCur[1][0]);

  double  *valuesDotDotCur1 = &(surf1->ValuesDotDotCur[0][0]);
  double  *valuesDotDotCur2 = &(surf1->ValuesDotDotCur[1][0]);


  double  af = surf1->td(2);
  double  d1 = surf1->td(5);
  double  aa = surf1->td(10);


   Idev2D(Idev);

   Klocal.setZero();
   Flocal.setZero();

   count = 1;   ll = 0;   err = 0;   isw = 3;
   dt = mpapTime.dt;

   count1 = 0;
   for(gp2=0;gp2<nGP2;gp2++)
   {
      for(gp1=0;gp1<nGP1;gp1++)
      {
          index = count1*2;

          //cout << " aaaaaaaaaaaa " << endl;
          //surf0->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, Jac);
          surf0->ShapeFunDerivatives2(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, d2N_dx2, d2N_dy2, d2N_dxy, d2N_dyx, Jac);
          //cout << " bbbbbbbbbbbbb " << endl;
          fact = gaussweights2[gp2] * gaussweights1[gp1] * thick * JacMultFact;

          dvol0 = Jac * fact;

          //cout << " bbbbbbbbbbbbb " << endl;
          surf1->deformationGradient(startindex[0], startindex[1], 1, dN_dx, dN_dy, F, detF);

          //surf1->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, Jac);
          surf1->ShapeFunDerivatives2(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, d2N_dx2, d2N_dy2, d2N_dxy, d2N_dyx, Jac);

          dvol = Jac * fact;

          veloCur[0] = veloCur[1] = 0.0;
          acceCur[0] = acceCur[1] = 0.0;

          pres = 0.0;
          dp[0] = dp[1] = 0.0;
          D2u.setZero();
          for(ii=0;ii<nlbf;ii++)
          {
            kk = tt[ii];

            bb1 = valuesCur1[kk];
            bb2 = valuesCur2[kk];
            bb3 = valuesCur3[kk];

            pres  += N[ii]*bb3;
            dp[0] += dN_dx[ii]*bb3;
            dp[1] += dN_dy[ii]*bb3;

            D2u(0,0) += bb1*d2N_dx2[ii];
            D2u(0,1) += bb1*d2N_dxy[ii];
            D2u(0,2) += bb1*d2N_dy2[ii];

            D2u(1,0) += bb2*d2N_dx2[ii];
            D2u(1,1) += bb2*d2N_dxy[ii];
            D2u(1,2) += bb2*d2N_dy2[ii];

            veloCur[0] += valuesDotCur1[kk]*N[ii];
            veloCur[1] += valuesDotCur2[kk]*N[ii];

            acceCur[0] += valuesDotDotCur1[kk]*N[ii];
            acceCur[1] += valuesDotDotCur2[kk]*N[ii];
          }

          F33 = 1.0;

          matlib2d_(matDat, F, &F33, stre, D11[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);
          if(err !=0)           return 1;

          dvol *= F33;

        //printf(" FFFFFFFF ");        printf("\t%12.8f\t%12.8f\t%12.8f\t%12.8f\t%12.8f\t%12.8f\n", F[0], F[1], F[2], F[3], detF, dvol);
        //printf(" stresses ");        printf("\t%12.8f\t%12.8f\t%12.8f\t%12.8f\n", stre[0], stre[1], stre[2], pres);

          pbar = (stre[0]+stre[1]+stre[2])/3.0;

          strdev[0] = stre[0] - pbar;
          strdev[1] = stre[1] - pbar;
          strdev[2] = stre[2] - pbar;
          strdev[3] = stre[3];

          stre[0] = strdev[0] + pres;
          stre[1] = strdev[1] + pres;
          stre[2] = strdev[2] + pres;

        //printf(" stresses ");        printf("\t%12.8f\t%12.8f\t%12.8f\t%12.8f\n\n", stre[0], stre[1], stre[2], pres);
//        printf(" volumes   ");        printf("\t%12.8f\t%12.8f\t%12.8f\t%12.8f\n\n", F33, detF, dvol0, dvol);
//          printf(" detF ");        printf("\t%12.8f\n", detF);

//          printf(" detF ");        printf("\t%12.8f\t%12.8f\t%12.8f\t%12.8f\n", detF,dvol0,dvol, dvol0*detF);

          for(ii=0;ii<4;ii++)
          {
             for(jj=0;jj<4;jj++)
             {
                cctmp[ii][jj] = 0.0;
                for(mm=0;mm<4;mm++)
                  cctmp[ii][jj] += Idev[ii][mm] * D11[mm][jj];
             }
          }

          for(ii=0;ii<4;ii++)
          {
             for(jj=0;jj<4;jj++)
             {
                D11[ii][jj] = 0.0;
                for(mm=0;mm<4;mm++)
                  D11[ii][jj] += cctmp[ii][mm] * Idev[mm][jj];
             }
          }


          fact = 2.0 * (pbar - pres);
          fact1 = (r2d3*pbar - pres);

          for(ii=0;ii<3;ii++)
          {
             for(jj=0;jj<4;jj++)
             {
                D11[ii][jj] -= r2d3*strdev[jj] ;
                D11[jj][ii] -= r2d3*strdev[jj] ;
             }
             for(jj=0;jj<3;jj++)
                D11[ii][jj] -= fact1;

             D11[ii][ii] += fact;
          }

          D11[3][3] += 0.5*fact;

          //==============================================
          // CALCULATE STIFFNESS matrices and RESIDUALs
          //==============================================

          dp[0] += 0.0*2.0*mu*(2.0*D2u(0,0)/3.0 + D2u(0,2)/2.0 + D2u(1,1)/6.0);
          dp[1] += 0.0*2.0*mu*(D2u(0,1)/6.0 + 2.0*D2u(1,2)/3.0 + D2u(1,0)/2.0);

          //force[0] = bforce[0] * rho0 ;
          //force[1] = bforce[1] * rho0 ;
          force[0] = 0.0;
          force[1] = 0.0;
          force[0] += rho0*acceCur[0];
          force[1] += rho0*acceCur[1];

          rStab[0] = rho0*acceCur[0] - dp[0];
          rStab[1] = rho0*acceCur[1] - dp[1];

          rStab[0] = - dp[0];
          rStab[1] = - dp[1];

          volstr = detF - 1.0 - pres*eps;

          for(ii=0;ii<nlbf;ii++)
          {
             TI   = 3*ii;
             TIp1 = TI+1;
             TIp2 = TI+2;

             bb1 = dN_dx[ii] * dvol;
             bb2 = dN_dy[ii] * dvol;
             bb3 = N[ii] * dvol;
             bb4 = N[ii] * dvol0;

             bc[0][0] = bb1 * D11[0][0] + bb2 * D11[3][0];
             bc[0][1] = bb1 * D11[0][1] + bb2 * D11[3][1];
             bc[0][2] = bb1 * D11[0][3] + bb2 * D11[3][3];

             bc[1][0] = bb2 * D11[1][0] + bb1 * D11[3][0];
             bc[1][1] = bb2 * D11[1][1] + bb1 * D11[3][1];
             bc[1][2] = bb2 * D11[1][3] + bb1 * D11[3][3];

             sig[0] = bb1 * stre[0] + bb2 * stre[3] ;
             sig[1] = bb1 * stre[3] + bb2 * stre[1] ;

             Flocal[TI]   -= (bb4*force[0] + sig[0]) ;
             Flocal[TIp1] -= (bb4*force[1] + sig[1]) ;
             Flocal[TIp2] -=  bb4*volstr;

             // PSPG stabilization
             Flocal[TIp2] -= tau*(bb1*rStab[0] + bb2*rStab[1]);

             acceFact = d1*rho0*bb4;

             for(jj=0;jj<nlbf;jj++)
             {
                TJ   = 3*jj;
                TJp1 = TJ+1;
                TJp2 = TJ+2;

                cc1 = dN_dx[jj];
                cc2 = dN_dy[jj];
                cc3 = N[jj];

                fact2 = acceFact*N[jj];

                Klocal(TI,   TJ)    += fact2 ;
                Klocal(TIp1, TJp1)  += fact2 ;

                fact = sig[0] * cc1 + sig[1] * cc2;

                Klocal(TI, TJ)      +=  af*(bc[0][0] * cc1 + bc[0][2] * cc2 + fact) ;
                Klocal(TI, TJp1)    +=  af*(bc[0][1] * cc2 + bc[0][2] * cc1) ;
                Klocal(TI, TJp2)    +=  af*(bb1 * cc3);

                Klocal(TIp1, TJ)    +=  af*(bc[1][0] * cc1 + bc[1][2] * cc2) ;
                Klocal(TIp1, TJp1)  +=  af*(bc[1][1] * cc2 + bc[1][2] * cc1 + fact) ;
                Klocal(TIp1, TJp2)  +=  af*(bb2 * cc3);

                Klocal(TIp2, TJ)    +=  af*(bb3 * cc1);
                Klocal(TIp2, TJp1)  +=  af*(bb3 * cc2);
                Klocal(TIp2, TJp2)  -=  af*(bb4 * cc3)*eps;

                // PSPG stabilization

                //Kpu(ii, twoJ)   -= tau*2.0*mu*( bb1 * (2.0*d2N_dx2[jj]/3.0+d2N_dy2[jj]/2.0) + bb2*(d2N_dxy[jj]/6.0));
                //Kpu(ii, twoJp1) -= tau*2.0*mu*( bb1 * (d2N_dxy[jj]/6.0) + bb2*(d2N_dx2[jj]/2.0+2.0*d2N_dy2[jj]/3.0));

                //dp[0] += 2.0*mu*(2.0*D2u(0,0)/3.0 + D2u(0,2)/2.0 + D2u(1,1)/6.0);
                //dp[1] += 2.0*mu*(D2u(0,1)/6.0 + 2.0*D2u(1,2)/3.0 + D2u(1,0)/2.0);

                //fact2 -= ( muTaf*d2N(jj) );
                fact2 = 0.0;
                //fact2 = d1*rho0*cc3;

                Dj(0,0) = fact2 - af*2.0*mu*(2.0*d2N_dx2[jj]/3.0+d2N_dy2[jj]/2.0)*0.0;
                Dj(0,1) = -af*2.0*mu*(d2N_dxy[jj]/6.0)*0.0;
                Dj(0,2) = -af*dN_dx[jj];
                Dj(1,0) = -af*2.0*mu*(d2N_dxy[jj]/6.0)*0.0;
                Dj(1,1) = fact2 - af*2.0*mu*(d2N_dx2[jj]/2.0+2.0*d2N_dy2[jj]/3.0)*0.0;
                Dj(1,2) = -af*dN_dy[jj];

                Klocal(TIp2, TJ)   += (bb1*Dj(0,0) + bb2*Dj(1,0))*tau;
                Klocal(TIp2, TJp1) += (bb1*Dj(0,1) + bb2*Dj(1,1))*tau;
                Klocal(TIp2, TJp2) += (bb1*Dj(0,2) + bb2*Dj(1,2))*tau;
             }
          }

          count++;
          count1++;
          ll += nivGP;
      }
   }

  return 0;
}
//




int NurbsElem2DStructMixed2fieldStabilised::calcStiffnessAndResidualAxsySS()
{
  return 0;
}



int NurbsElem2DStructMixed2fieldStabilised::calcStiffnessAndResidualAxsyFS()
{
  return 0;
}


int NurbsElem2DStructMixed2fieldStabilised::calcInternalForces()
{

  return 0;
}




int NurbsElem2DStructMixed2fieldStabilised::calcOutput(double u1, double v1)
{

  return 0;
}







void NurbsElem2DStructMixed2fieldStabilised::discreteContourplot(int vartype, int varindex, int index, int nCol, double umin, double umax)
{
  if(index > nivGP)
  {
     cout << '\t' << " Error in NurbsElem2DStructMixed2fieldStabilised::contourplot " << endl;
     return;
  }

  if(varindex > 4)
    varindex -= 2;


   double outval[500];

   switch(vartype)
   {
       case 0:  // plot total strain
       case 1:  // plot elastic strain
       case 2:  // plot plastic strain

                projectStrain(vartype, varindex, outval);

              break;

       case 3:  // plot stress

                projectStress(varindex, outval);

              break;

       case 4:  // plot element internal variables

                projectIntVar(index, outval);

              break;

       default:

              cout  << "           Invalid Variable Type to project " << endl;
              break;

   }

  double uu, vv, du, dv;

  du = uvalues[2]/nGP1;
  dv = vvalues[2]/nGP2;

  ListArray<EPOINT> S1;
  S1.setDim( (nGP1+1)*(nGP2+1) );

  int count=0, ii, jj;

  vv = vvalues[0];

  if(finite)
  {
     for(jj=0;jj<=nGP2;jj++)
     {
        uu = uvalues[0];
        for(ii=0;ii<=nGP1;ii++)
        {
           S1[count] = surf1->SurfacePoint(uu, vv).CalcEuclid();
           count++;
           uu += du;
        }
        vv += dv;
     }
  }
  else
  {
     for(jj=0;jj<=nGP2;jj++)
     {
        uu = uvalues[0];
        for(ii=0;ii<=nGP1;ii++)
        {
           S1[count] = surf0->SurfacePoint(uu, vv).CalcEuclid();
           count++;
           uu += du;
        }
        vv += dv;
     }

  }

  EPOINT *EP;

  double   x1[2], x2[2], x3[2], x4[2], u1;

  int ind1, ind2, nGP1p1;

  count=0;
  nGP1p1 = nGP1+1;
  for(jj=0;jj<nGP2;jj++)
  {
      ind1 = (nGP1p1)*jj;
      ind2 = (nGP1p1)*(jj+1);

      for(ii=0;ii<nGP1;ii++)
      {
          u1 = outval[count];

          EP = &(S1[ind1+ii]);
          x1[0] = EP->x; x1[1] = EP->y;

          EP = &(S1[ind2+ii]);
          x2[0] = EP->x; x2[1] = EP->y;

          EP = &(S1[ind2+ii+1]);
          x3[0] = EP->x; x3[1] = EP->y;

          EP = &(S1[ind1+ii+1]);
          x4[0] = EP->x; x4[1] = EP->y;

          // contour plot for 1st triangle
          //plot.triangleContourPlot(x1, x2, x3, u1, u1, u1, umin, umax, nCol);

          // contour plot for 2nd triangle
          //plot.triangleContourPlot(x1, x3, x4, u1, u1, u1, umin, umax, nCol);

          count++;
      }
  }

  return;
}








void NurbsElem2DStructMixed2fieldStabilised::projectToKnots(bool extrapolateFlag, int vartype, int varindex, int index)
{
/*
   vals2project[0] = intVar2[indx];
   vals2project[1] = intVar2[(nGP1-1)*nivGP+indx];
   vals2project[2] = intVar2[nGP1*(nGP2-1)*nivGP+indx];
   vals2project[3] = intVar2[(nGP1*nGP2-1)*nivGP+indx];
*/

   double outval[500];

   switch(vartype)
   {
       case 1:  // plot total strain
       case 2:  // plot elastic strain
       case 3:  // plot plastic strain

                projectStrain(vartype, varindex, outval);

              break;

       case 4:  // plot stress

                projectStress(varindex, outval);

              break;

       case 5:  // plot element internal variables

                projectIntVar(index, outval);

              break;

       default:

              cout  << "           Invalid Variable Type to project " << endl;
              break;

   }


   assert(vals2project.n == 4);

    if(extrapolateFlag)
    {
       for(int ii=0;ii<4;ii++)
         vals2project[ii] = extrapolate(nGP1, (ii+1), outval);
    }
    else
    {
       vals2project[0] = outval[0];
       vals2project[1] = outval[nGP1-1];
       vals2project[2] = outval[nGP1*(nGP2-1)];
       vals2project[3] = outval[nGP1*nGP2-1];
    }

   //cout << '\t' << vals2project << endl; cout << endl;

  return;
}




void NurbsElem2DStructMixed2fieldStabilised::projectStress(int varindex, double* outval)
{
//
   int nivEL = nGP * nivGP;
   for(int ii=0;ii<nivEL;ii++)
     intVar2[ii] = intVar1[ii];

   int  err, isw, count, count1, ll, gp1, gp2, index;

   double  F[4], detF, F33, fact, dt, Jac, pres;

   double  cc[4][4], stre[4], N[nlbf], dN_dx[nlbf], dN_dy[nlbf];

   count  = 1;   ll = 0;  err = 0;   isw  = 3;
   dt = mpapTime.dt;

   count1 = 0;
   for(gp2=0;gp2<nGP2;gp2++)
   {
      for(gp1=0;gp1<nGP1;gp1++)
      {
          index = count1*2;

          surf0->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, Jac);

          surf1->deformationGradient(startindex[0], startindex[1], 1, dN_dx, dN_dy, F, detF);

          F33 = 1.0;

          matlib2d_(matDat, F, &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);

          pres = surf1->computeValue(3, knotsAtGPs[index], knotsAtGPs[index+1]);

          //cout << " pres " << pres << endl;

          fact = pres - (stre[0]+stre[1]+stre[2])/3.0 ;

          stre[0] += fact;
          stre[1] += fact;
          stre[2] += fact;

          if(varindex < 4)
             outval[count1] = stre[varindex];
          else if(varindex == 4)
             outval[count1] = sqrt((pow(stre[0]-stre[1],2.0) + pow(stre[1]-stre[2], 2.0) + pow(stre[2]-stre[0], 2.0) + 6.0* stre[3]*stre[3])/2.0);
          else if(varindex == 5)
             outval[count1] = pres;
          else
          {
             cout << '\t' << "    NurbsElem2DStructMixed2fieldStabilised::projectStress .... : Error in 'varindex' " << endl;
             return;
          }

          count++;
          count1++;
          ll += nivGP;
     }
   }
//
  return;
}




void NurbsElem2DStructMixed2fieldStabilised::projectStrain(int vartype, int varindex, double* outval)
{
/*
   int nivEL = nGP * nivGP;
   for(int ii=0;ii<nivEL;ii++)
     intVar2[ii] = intVar1[ii];

   bool  intVarFlag = (nivGP > 0);



   // loop over Gauss points
   for(gp2=0;gp2<nGP2;gp2++)
   {
      for(gp1=0;gp1<nGP1;gp1++)
      {

        switch(vartype)
        {
             case 1: // total strain components

                     switch(varindex)
                     {
                          case 0:    // epsxx

                              outval[count1] = F[0] - 1.0;           break;

                          case 1:    // epsyy

                              outval[count1] = F[3] - 1.0;           break;

                          case 2:    // epsxy

                              outval[count1] = F[1] + F[2];           break;
                     }

               break;

             case 2: // plastic strain components

                 if(intVarFlag)
                 {
                     switch(varindex)
                     {
                          case 0:    // epspxx

                              outval[count1] = intVar2[0];         break;

                          case 1:    // epspyy

                              outval[count1] = intVar2[1];         break;

                          case 2:    // epspxy

                              outval[count1] = intVar2[3];         break;

                          case 3:    // epspeqv

                              outval[count1] = intVar2[6];         break;
                     }
                 }

               break;

             case 3: // elastic strain components


               if(intVarFlag)
               {
                     switch(varindex)
                     {
                          case 0:    // epspxx

                              outval[count1] = F[0] - 1.0 - intVar2[0];         break;

                          case 1:    // epspyy

                              outval[count1] = F[3] - 1.0 - intVar2[1];         break;

                          case 2:    // epspxy

                              outval[count1] = F[1] + F[2] - intVar2[3];         break;
                     }

               }
               else
               {
                     switch(varindex)
                     {
                          case 0:    // epsxx

                              outval[count1] = F[0] - 1.0;           break;

                          case 1:    // epsyy

                              outval[count1] = F[3] - 1.0;           break;

                          case 2:    // epsxy

                              outval[count1] = F[1] + F[2];           break;
                     }
               }

               break;
        }

        count++;
        count1++;
        ll += nivGP;

    }
  }

*/

return;
}






void NurbsElem2DStructMixed2fieldStabilised::projectIntVar(int index, double* outval)
{
   int ind1, ind2, ii, jj;

   ind1 = 0;
   for(jj=0;jj<nGP2;jj++)
   {
       for(ii=0;ii<nGP1;ii++)
       {
           outval[ind1] = intVar2[ind1*nivGP+index];
           ind1++;
       }
   }


   return;
}


/*
void NurbsElem2DStructMixed2fieldStabilised::AssembleElementMatrix2(int index, MatrixXd& stiff, MatrixXd& CC)
{
    int row, aa, bb, *colC,  ind;

    ind  = surf2->IEN[elenum].n;
    colC = &(surf2->IEN[elenum][0]);

    for(aa=0;aa<nsize;aa++)
    {
        row = forassy[aa];
        for(bb=0;bb<nsize;bb++)
            stiff(row, forassy[bb]) += stiffness_local[aa][bb];

        for(bb=0;bb<ind;bb++)
           CC(row, colC[bb]) += Kup(aa,bb);
    }

  return;
}
*/

void NurbsElem2DStructMixed2fieldStabilised::AssembleElementMatrix2(int index, MatrixXd& mat, MatrixXd& CC)
{
    int *tt1, *tt2;
    int row, col, aa, bb, ind, n1, n2;
    
    //cout <<  surf0->LM[elenum] << endl;
    //cout << " \n\n " << endl;
    //cout <<  surf2->LM[elenum] << endl;
    //cout << " \n\n " << endl;


    tt1 = &(surf0->LM[elenum][0]);
    n1  = surf0->LM[elenum].n;
    tt2 = &(surf2->LM[elenum2][0]);
    n2  = surf2->LM[elenum2].n;

    if(index == 1)
    {
      for(aa=0;aa<n1;aa++)
      {
        row = tt1[aa];
        if(row != -1)
        {
          for(bb=0;bb<n1;bb++)
          {
            col = tt1[bb];
            if(col != -1)
              mat(row, col) += Klocal(aa, bb);
          }
        }
      }
    }
    if(index == 2)
    {
      for(aa=0;aa<n1;aa++)
      {
        row = tt1[aa];
        if(row != -1)
        {
          for(bb=0;bb<n2;bb++)
          {
            col = tt2[bb];
            if(col != -1)
              mat(row, col) += Kup(aa,bb);
          }
        }
      }
    }
    if(index == 3)
    {
      for(aa=0;aa<n2;aa++)
      {
        row = tt2[aa];
        if(row != -1)
        {
          for(bb=0;bb<n2;bb++)
          {
            col = tt2[bb];
            if(col != -1)
              mat(row, col) += Kpp(aa,bb);
          }
        }
      }
    }

    return;
}


void NurbsElem2DStructMixed2fieldStabilised::AssembleElementVector2(bool firstIter, int flag, VectorXd& rhs, double* reac, VectorXd& uu)
{
   // flag == true  ---> just external force vector
   // flag == false ---> internal load vector

   if(flag == 1)
   {
      for(int aa=0;aa<nsize;aa++)
      {
         rhs(forassy[aa]) += resi[aa];
      }
   }
   else if(flag == 0)
   {
      for(int aa=0;aa<nsize;aa++)
      {
         rhs(forassy[aa])  += resi[aa];
         reac[forassy[aa]] += resi[aa];
      }
   }
   else if(flag == 2)
   {
      int  size2 = surf2->nsize;
      int  *tt2  = &(surf2->LM[elenum][0]);
      int aa, bb, ind;

      for(aa=0;aa<size2;aa++)
      {
          ind = tt2[aa];

          for(bb=0;bb<nsize;bb++)
          {
             rhs[ind] += Kup(bb,aa) * uu(forassy[bb]);
          }
      }
   }
   else
   {
      cerr << " ERROR in flag in 'NurbsElem2DStructMixed2fieldStabilised::AssembleElementVector2' " << endl;
      return;
   }


//  cout << " resi " << resi << endl;

  return;
}






void NurbsElem2DStructMixed2fieldStabilised::AssembleElementMatrix(int index, MatrixSparseArray<double>& mtx)
{
    int nn=0, aa, bb, size2;

    size2 = surf2->nsize;

    //cout << " hhhhhhhhh " << endl;
    for(aa=0;aa<nsize;aa++)
    {
      for(bb=0;bb<nsize;bb++)
      {
        nn = forassembly[aa][bb];
        if(nn != -1)
          mtx.x[nn-1] += Klocal(aa, bb);
      }

      for(bb=0;bb<size2;bb++)
      {
        nn = forassyKup[aa][bb];
        if(nn != -1)
          mtx.x[nn-1] += Kup(aa,bb);

        nn = forassyKpu[bb][aa];
        if(nn != -1)
          mtx.x[nn-1] += Kup(aa,bb);
      }
    }
    for(aa=0;aa<size2;aa++)
    {
      for(bb=0;bb<size2;bb++)
      {
        nn = forassyKtt[aa][bb];
        if(nn != -1)
          mtx.x[nn-1] += Kpp(aa,bb);
      }
    }

  return;
}





void NurbsElem2DStructMixed2fieldStabilised::AssembleElementVector(bool firstIter, bool flag, double* rhs, double* reac, int start1, int start2)
{
   // flag == true  ---> just external force vector
   // flag == false ---> internal load vector + contributions from nodes with specified displacement BCs

   int *tt;

   //cout << " primvar " << endl;
   //cout << primvar << endl;

    tt = &(surf0->LM[elenum][0]);

    if(flag)
    {
      for(int aa=0;aa<nsize;aa++)
      {
         if(tt[aa] != -1)
            rhs[tt[aa]] += Flocal[aa];
      }
    }
    else
    {
      double fact;
      int aa, bb, ind;

      //cout << " kkkkkkkkkkkk " << endl;

      for(aa=0;aa<nsize;aa++)
      {
         if(tt[aa] != -1)
           rhs[tt[aa]] += Flocal[aa];

         // add up reaction forces
         reac[forassy[aa]] += Flocal[aa];
      }
      //cout << " qqqqqqqqqqqq " << endl;

      // contribution to the pressure variables from the applied displacements

      if(firstIter)
      {
        for(aa=0;aa<nsize;aa++)
        {
            if(tt[aa] == -1)
            {
                fact = mpapTime.dt * primvar[aa];
                for(bb=0;bb<nsize;bb++)
                {
                    if(tt[bb] != -1)
                    {
                      //printf("\t%5d\t%5d\t%5d\t%12.8f\t%12.8f\n\n", aa, bb, tt[bb], stiffness_local[bb][aa], fact);
                      rhs[tt[bb]] -= Klocal(bb, aa) * fact;
                    }
                }
            }
        }
      }
   }

//  cout << " resi " << resi << endl;

  tt= NULL;

  return;
}


/*
void  NurbsElem2DStructMixed2fieldStabilised::AssembleElementMatrix(int index, Mat mtx, int start1, int start2)
{
    PetscErrorCode ierr;
    int  ii, jj, nn=0, aa, bb, size2, ind, *tt1, *tt2;

    tt1 = &(surf0->LM[elenum][0]);
    tt2 = &(surf2->LM[elenum2][0]);
    size2 = surf2->nsize;
    
    //cout << elenum << '\t' << elenum2 << endl;
    //cout << nsize << '\t' << size2 << endl;
    //cout << surf0->LM[elenum] << endl;
    //cout << surf2->LM[elenum2] << endl;

    for(ii=0;ii<nsize;ii++)
    {
      aa = tt1[ii];
      if( aa != -1)
      {
        for(jj=0;jj<nsize;jj++)
        {
          bb = tt1[jj];
          if(bb != -1)
          {
            //cout << ii << '\t' << jj << '\t' << aa << '\t' << bb << endl;
            ierr = MatSetValues(mtx, 1, &aa, 1, &bb, &(Klocal(ii, jj)), ADD_VALUES);
          }
        }
        //cout << " ii = " << ii << endl;

        for(jj=0;jj<size2;jj++)
        {
          bb = tt2[jj];
          if(bb != -1)
          {
            bb += start1;
            ierr = MatSetValues(mtx, 1, &aa, 1, &bb, &(Kup(ii, jj)), ADD_VALUES);
            ierr = MatSetValues(mtx, 1, &bb, 1, &aa, &(Kup(ii, jj)), ADD_VALUES);
          }
        }
        //cout << " ii = " << ii << endl;
      }
    }
    //cout << " qqqqqqqqqqqq " << endl;
    for(ii=0;ii<size2;ii++)
    {
      aa = tt2[ii];
      if( aa != -1)
      {
        aa += start1;
        for(jj=0;jj<size2;jj++)
        {
          bb = tt2[jj];
          if(bb != -1)
          {
            bb += start1;
            ierr = MatSetValues(mtx, 1, &aa, 1, &bb, &(Kpp(ii, jj)), ADD_VALUES);
          }
        }
      }
    }

  return;
}
*/


void  NurbsElem2DStructMixed2fieldStabilised::AssembleElementMatrix(int index, SparseMatrixXd& mtx, int start1, int start2)
{
    int  ii, jj, nn=0, aa, bb, size2, ind, *tt1, *tt2;

    tt1 = &(surf0->LM[elenum][0]);
    
    /*
    cout << elenum << '\t' << elenum2 << endl;
    cout << nsize << '\t' << size2 << endl;
    cout << surf0->LM[elenum] << endl;

   printMatrix(Klocal);
   printf("\n\n\n");
   */


    for(ii=0;ii<nsize;ii++)
    {
      aa = tt1[ii];
      if( aa != -1)
      {
        for(jj=0;jj<nsize;jj++)
        {
          bb = tt1[jj];
          if(bb != -1)
          {
            //cout << ii << '\t' << jj << '\t' << aa << '\t' << bb << endl;
            mtx.coeffRef(aa, bb) += Klocal(ii, jj);
          }
        }
        //cout << " ii = " << ii << '\t' << aa << endl;
      }
    }

  return;
}



void NurbsElem2DStructMixed2fieldStabilised::toPostprocess(int vartype, int varindex, int type, SparseMatrixXd&  coeffMat, VectorXd& rhsVec)
{
/*
   double F[4], detF=0.0, F33, Jac, dt, dN_dx[nlbf], dN_dy[nlbf], stre[4], cc[4][4];

   int   err,  isw,  count,  count1, index, ll = 0, ii, jj, gp1, gp2, row, col;

   MatrixXd  Nlocal(nlbf,nlbf);
   VectorXd  NN(nlbf), rhslocal(nlbf);

   Nlocal.setZero();   
   rhslocal.setZero();

   double *gaussweights = &(surf0->gaussweights[0]);

   count = 1;   ll = 0;   err = 0;   isw = 3;
   dt = mpapTime.dt;

   count1 = 0;
   for(gp2=0;gp2<nGP2;gp2++)
   {
   for(gp1=0;gp1<nGP1;gp1++)
   {
        index = 2*count1;

        surf0->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), &NN(0), dN_dx, dN_dy, Jac);

        surf1->deformationGradient(startindex[0], startindex[1], 1, dN_dx, dN_dy, F, detF);

        if(sss == 1)  // plane stress
        {
          if(finite)
            F33 = 1.0/sqrt(detF);
          else
            F33 = 3.0 - F[0] - F[3];
        }
        else if(sss == 2)    // plane strain
          F33 = 1.0;

        matlib2d_(matDat, F, &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);
        count++;
        count1++;
        ll += nivGP;
        
        Nlocal += NN * NN.transpose();
        
        if(varindex < 4)
           rhslocal += NN * stre[varindex];
        else if(varindex == 4)
           rhslocal += NN * sqrt((pow(stre[0]-stre[1],2.0) + pow(stre[1]-stre[2], 2.0) + pow(stre[2]-stre[0], 2.0) + 6.0 * stre[3]*stre[3])/2.0);
        else if(varindex == 5)
           rhslocal += NN * (stre[0]+stre[1]+stre[2])/3.0;

  }//gp1
  }//gp2


    int *tt;
    tt = &(surf0->IEN[elenum][0]);
      
    for(ii=0;ii<nlbf;ii++)
    {
       row = tt[ii];
       rhsVec(row) += rhslocal(ii);
       for(jj=0;jj<nlbf;jj++)
       {
          col = tt[jj];
          coeffMat.coeffRef(row, col) += Nlocal(ii, jj);
       }
    }
*/

  return;
}







int NurbsElem2DStructMixed2fieldStabilised::calcError(int ind)
{
   int ii, jj, gp1, gp2, count1,  TI, index;
   
   double rad, theta, Jac, fact, dvol0;
   
   double  N[nlbf], dN_dx[nlbf], dN_dy[nlbf], v1, v2, disp1[2], disp2[2];
   
   //PlateWithHole  analy(sss, matDat[1], matDat[2]);
   ThickCylinder  analy(sss, matDat[1], matDat[2]);
   //ElasticityEx2  analy(matDat[0], matDat[1]);

   EPOINT  EP, EP2;

   double *gaussweights1 = &(surf0->gaussweights1[0]);
   double *gaussweights2 = &(surf0->gaussweights2[0]);

   double  *values1 = &(surf1->Values[0][0]);
   double  *values2 = &(surf1->Values[1][0]);
   double  *values3 = &(surf1->Values[2][0]);

   int *tt = &(surf0->IEN[elenum][0]);

   elemError = 0.0;

   if(ind == 1) // x - displacement
   {
      count1 = 0;
      for(gp2=0;gp2<nGP2;gp2++)
      {
      for(gp1=0;gp1<nGP1;gp1++)
      {
          index = count1*2;

          surf0->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, Jac);
      
          fact = gaussweights2[gp2] * gaussweights1[gp1] * thick * JacMultFact;

          dvol0 = Jac * fact;
          count1++;

          EP  = surf0->SurfacePoint(knotsAtGPs[index], knotsAtGPs[index+1]).CalcEuclid();

          EP2 = surf1->SurfacePoint(knotsAtGPs[index], knotsAtGPs[index+1]).CalcEuclid();

          rad = sqrt(EP.x*EP.x + EP.y*EP.y) ;

          theta = atan2(EP.y, EP.x);

          disp1[0] = analy.dispX(rad, theta);
          disp1[1] = analy.dispY(rad, theta);

          //disp1[0] = analy.dispX(EP.x, EP.y);
          //disp1[1] = analy.dispY(EP.x, EP.y);

          disp2[0] = disp2[1] = 0.0;
          //for(ii=0;ii<nlbf;ii++)
          //{
            //disp2[0] += values1[tt[ii]]* N[ii];
            //disp2[1] += values2[tt[ii]]* N[ii];
          //}
          disp2[0] = EP2.x - EP.x;
          disp2[1] = EP2.y - EP.y;
          
          //v2 = disp[0]*cos(theta) ;

          //printf(" \t %12.8f \t %12.8f \n ", EP.x, EP.y);
          //printf(" \t %12.8f \t %12.8f \t %12.8f \n ", v1, v2, (-disp[0]*sin(theta) + disp[1]*cos(theta)));

          //v1 -= v2;
          disp1[0] -= disp2[0] ;
          disp1[1] -= disp2[1] ;

          fact = disp1[0]*disp1[0] + disp1[1]*disp1[1];
          //fact = er*er ;
          
          elemError += ( fact * dvol0 );
      }
      }
   }
   if(ind == 2) // y - displacement
   {

   }
   if(ind == 3) // stress
   {

      int ll, isw, err, count;

      double F[4], detF, F33, stre[4], stre2[4], cc[4][4], bb1, bb2, dt;
      int  sizep = surf2->nlbf;
      double  pres, dummy, volstr;

      count = 1;   ll = 0;   err = 0;   isw = 3;
      dt = mpapTime.dt;

      count1 = 0;
      for(gp2=0;gp2<nGP2;gp2++)
      {
      for(gp1=0;gp1<nGP1;gp1++)
      {
          index = count1*2;

          surf0->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, Jac);

          fact = gaussweights2[gp2] * gaussweights1[gp1] * thick * JacMultFact;

          dvol0 = Jac * fact;

          EP = surf0->SurfacePoint(knotsAtGPs[index], knotsAtGPs[index+1]).CalcEuclid();

          rad = sqrt(EP.x*EP.x + EP.y*EP.y) ;

          theta = atan2(EP.y, EP.x);

          stre[0] = analy.stressXX(rad, theta);
          stre[1] = analy.stressYY(rad, theta);
          stre[2] = 0.0;
          if(sss == 2)
            stre[2] = matDat[2]*(stre[0]+stre[1]);
          stre[3] = analy.stressXY(rad, theta);

          F[0] = F[1] = F[2] = F[3] = 0.0;
          pres = 0.0;
          for(ii=0;ii<nlbf;ii++)
          {
            jj = tt[ii];

            bb1 = values1[jj];
            bb2 = values2[jj];

            F[0] += bb1*dN_dx[ii];
            F[2] += bb1*dN_dy[ii];
            F[1] += bb2*dN_dx[ii];
            F[3] += bb2*dN_dy[ii];

            pres += values3[ii] * N[ii];
          }
          F[0] += 1.0;
          F[3] += 1.0;
          
          //cout << " pres = " << pres << endl;
          
          //surf1->deformationGradient(startindex[0], startindex[1], 1, dN_dx, dN_dy, F, detF);

          volstr = (F[0]+F[3]-2.0);

          F33 = 1.0;

          matlib2d_(matDat, F, &F33, stre2, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);

          if(err !=0)           return 1;

          count++;
          count1++;
          ll += nivGP;

          dummy = pres - (stre2[0]+stre2[1]+stre2[2])/3.0 ;

          stre2[0] += dummy;
          stre2[1] += dummy;
          stre2[2] += dummy;

          for(ii=0;ii<4; ii++)
            stre[ii] -= stre2[ii];

        //printf("stresses... \t%20.18f\t%20.18f\t%20.18f\t%20.18f \n", stre[0], stre[1], stre[2], stre[3]);

          fact = stre[0]*stre[0] + stre[1]*stre[1] + stre[3]*stre[3] ;//+ stre[2]*stre[2] ;
          //fact = stre[0]*stre[0] ;

          elemError += ( fact * dvol0 );
      }
      }
   }
   if(ind == 4) // Energy norm
   {
      int ll, isw, err, count;

      double F1[4], F2[4], eps[4], detF, F33, stre[4], cc[4][4], bb1, bb2, dt;
      int  sizep = surf2->nlbf;
      double  pres, dummy, Nbar[sizep], volstr;

      count = 1;   ll = 0;   err = 0;   isw = 3;
      dt = mpapTime.dt;

      count1 = 0;
      for(gp2=0;gp2<nGP2;gp2++)
      {
      for(gp1=0;gp1<nGP1;gp1++)
      {
          index = count1*2;

          surf0->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, Jac);
      
          fact = gaussweights2[gp2] * gaussweights1[gp1] * thick * JacMultFact;

          dvol0 = Jac * fact;

          EP = surf0->SurfacePoint(knotsAtGPs[index], knotsAtGPs[index+1]).CalcEuclid();

          rad = sqrt(EP.x*EP.x + EP.y*EP.y) ;

          theta = atan2(EP.y, EP.x);

          F1[0] = analy.strainXX(rad, theta) ;
          F1[1] = analy.strainXY(rad, theta) ;
          F1[2] = F1[1];
          F1[3] = analy.strainYY(rad, theta) ;

          //stre[0] = analy.stressXX(EP.x, EP.y);
          //stre[1] = analy.stressYY(EP.x, EP.y);
          //stre[2] = 0.0;
          //stre[3] = analy.stressXY(EP.x, EP.y);

          //surf1->deformationGradient(startindex[0], startindex[1], 1, dN_dx, dN_dy, F, detF);
          
          //
          F2[0] = F2[1] = F2[2] = F2[3] = 0.0;
          pres = 0.0;

          for(ii=0;ii<nlbf;ii++)
          {
            jj = tt[ii];

            bb1 = values1[jj];
            bb2 = values2[jj];

            F2[0] += bb1*dN_dx[ii];
            F2[2] += bb1*dN_dy[ii];
            F2[1] += bb2*dN_dx[ii];
            F2[3] += bb2*dN_dy[ii];

            pres += values3[ii] * N[ii];
          }
          //F2[0] += 1.0;
          //F2[3] += 1.0;

        for(ii=0;ii<4; ii++)
          F1[ii] -= F2[ii];

        F1[0] += 1.0;
        F1[3] += 1.0;

        //printf("F... \t%20.18f\t%20.18f\t%20.18f\t%20.18f\t%20.18f \n", F[0], F[1], F[2], F[3], detF);

        // ADJUST F33 fOR 2D PROBLEMS BASED ON THE ASSUMPTIONS OF PLANE STRESS/PLANE STRAIN/AXISYMMETRIC

        if(sss == 1)  // plane stress
        {
          if(finite)
            F33 = 1.0/sqrt(detF);
          else
            F33 = 3.0 - F1[0] - F1[3];
        }
        else if(sss == 2)    // plane strain
          F33 = 1.0;


        matlib2d_(matDat, F1, &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);
        count++;
        count1++;
        ll += nivGP;

        dummy = pres - (stre[0]+stre[1]+stre[2])/3.0 ;

        stre[0] += dummy;
        stre[1] += dummy;
        stre[2] += dummy;

        eps[0] = F1[0]-1.0;
        eps[1] = F1[3]-1.0;
        eps[2] = 0.0;
        eps[3] = F1[1];

        fact = 0.5*(stre[0]*eps[0] + stre[1]*eps[1] + stre[3]*eps[3] );// + stre[3]*stre[3] ;

        elemError += ( fact * dvol0 );
      }
      }
   }
   if(ind == 5) // LS functional
   {

   }

  //if( (int) matDat[4] == 1)
    //printf(" \t element = %5d ... \t ... elemError  =   %12.6E \n " , elenum, elemError);

  return 0;
}


