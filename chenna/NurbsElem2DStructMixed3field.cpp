
#include <math.h>
#include "Debug.h"
//#include "Plot.h"
//#include "FunctionsElement.h"
#include "MpapTime.h"
#include "NurbsElem2DStructMixed3field.h"
#include "NurbsShapeFunctions.h"
#include <assert.h>
#include "ComputerTime.h"


using namespace std;

//extern Plot     plot;
extern MpapTime mpapTime;
extern ComputerTime       computerTime;



NurbsElem2DStructMixed3field::NurbsElem2DStructMixed3field()
{
  if (debug) cout << " constructor NurbsElem2DStructMixed3field\n\n";

//   cout << " constructor NurbsElem2DStructMixed3field\n\n";

  calcExtraMatrices = true;
}



NurbsElem2DStructMixed3field::~NurbsElem2DStructMixed3field()
{
  if (debug) cout << " destructor NurbsElem2DStructMixed3field\n\n";

  // cout << " destructor NurbsElem2DStructMixed3field\n\n";

}

int NurbsElem2DStructMixed3field::calcStiffnessMatrix(double dt)
{

  return 0;
}



void NurbsElem2DStructMixed3field::contourplot(int index, int nCol, double umin, double umax)
{
  return;
}

int NurbsElem2DStructMixed3field::calcMassMatrix(int lumpInd, double dt)
{


  return 0;
}



int NurbsElem2DStructMixed3field::calcStiffnessAndResidual()
{
   if(axsy)
   {
     cout << " AXISYMMETRIC case is not supported by NurbsElem2DStructMixed3field  " << endl;
     return -1;
   }

  if(finite)
    NurbsElem2DStructMixed3field::calcStiffnessAndResidual2();
  else
    NurbsElem2DStructMixed3field::calcStiffnessAndResidual1();

  return 0;
}



int NurbsElem2DStructMixed3field::calcStiffnessAndResidual1()
{
/*
   double  BULK = matDat[0];

   int  err, isw, count, count1, ll, ii, jj, twoI, twoJ, twoIp1, twoJp1, index, gp1, gp2, mm;

   int  sizeM = surf2->nlbf;

   double  F[4], detF, F33, dvol0, dvol, dt, Jac, pres, bb1, bb2, bb3, r1d3=1.0/3.0, r2d3=2.0/3.0;

   double  cc[4][4], stre[4], bc[2][3], N[nlbf], dN_dx[nlbf], dN_dy[nlbf], Nbar[sizeM];

   double fact, fact1, fact2, fact3, fact4, utemp, vtemp, volstr, volstr1, dummy, pbar;

   double  Idev[4][4], cctmp[4][4];

   Idev2D(Idev);

   if(calcExtraMatrices)
   {
      Kup.resize(nsize, sizeM);
      Ktt.resize(sizeM, sizeM);
      Ktp.resize(sizeM, sizeM);

      Kup.setZero();
      Ktt.setZero();
      Ktp.setZero();
   }


   double *gaussweights = &(surf0->gaussweights[0]);

   for(ii=0;ii<nsize;ii++)
     stiffness_local[ii].zero();

   resi.zero();
   resi2.zero();

   count = 1;   ll = 0;   err = 0;   isw = 3;
   dt = mpapTime.dt;

   count1 = 0;
   for(gp2=0;gp2<nGP2;gp2++)
   {
      for(gp1=0;gp1<nGP1;gp1++)
      {
          index = count1*2;

          utemp = knotsAtGPs[index];
          vtemp = knotsAtGPs[index+1];

          surf0->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, Jac);

          dvol0 = Jac * gaussweights[count1] * JacMultFact;

          surf1->deformationGradient(startindex[0], startindex[1], 1, dN_dx, dN_dy, F, detF);

          pres = surf2->computeValueAndShanpeFns(2, utemp, vtemp, Nbar);

          volstr1  = surf2->computeValue(1,utemp,vtemp);

          // modify F

          volstr = (F[0]+F[3]-2.0);

          fact = (volstr1 - volstr)/3.0;

          F[0] += fact;
          F[3] += fact;
          F33 = 1.0 + fact;

          matlib2d_(matDat, F, &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);
          if(err !=0)           return 1;

          pbar = (stre[0]+stre[1]+stre[2])/3.0;

          fact = pres - pbar;

          stre[0] += fact;
          stre[1] += fact;
          stre[2] += fact;

//        printf(" stresses ");        printf("\t%12.8f\t%12.8f\t%12.8f\t%12.8f\n\n", stre[0], stre[1], stre[2], pres);
//        printf(" values   ");        printf("\t%12.8f\t%12.8f\t%12.8f\t%12.8f\t%12.8f\n\n", detF, detFbar, pbar, phat, dvol);

          for(ii=0;ii<4;ii++)
          {
             for(jj=0;jj<4;jj++)
             {
                cctmp[ii][jj] = 0.0;
                for(mm=0;mm<4;mm++)
                  cctmp[ii][jj] += Idev[ii][mm] * cc[mm][jj];
             }
          }

          for(ii=0;ii<4;ii++)
          {
             for(jj=0;jj<4;jj++)
             {
                cc[ii][jj] = 0.0;
                for(mm=0;mm<4;mm++)
                  cc[ii][jj] += cctmp[ii][mm] * Idev[mm][jj];
             }
          }

          for(ii=0;ii<4;ii++)
          {
             stre[ii] *= dvol0;
             for(jj=0;jj<4;jj++)
               cc[ii][jj] *= dvol0;
          }

          //==============================================
          // CALCULATE TANGENT STIFFNESS and RESIDUAL
          //==============================================

          for(ii=0;ii<nlbf;ii++)
          {
             bb1 = dN_dx[ii];
             bb2 = dN_dy[ii];

             bc[0][0] = (bb1 * cc[0][0] + bb2 * cc[3][0]);
             bc[0][1] = (bb1 * cc[0][1] + bb2 * cc[3][1]);
             bc[0][2] = (bb1 * cc[0][3] + bb2 * cc[3][3]);

             bc[1][0] = (bb2 * cc[1][0] + bb1 * cc[3][0]);
             bc[1][1] = (bb2 * cc[1][1] + bb1 * cc[3][1]);
             bc[1][2] = (bb2 * cc[1][3] + bb1 * cc[3][3]);

             twoI   = 2*ii;
             twoIp1 = twoI+1;

             resi[twoI]   -= (bb1 * stre[0] + bb2 * stre[3]);
             resi[twoIp1] -= (bb1 * stre[3] + bb2 * stre[1]) ;

             for(jj=ii;jj<nlbf;jj++)
             {
                twoJ   = 2*jj;
                twoJp1 = twoJ+1;

                bb1 = dN_dx[jj];
                bb2 = dN_dy[jj];

                stiffness_local[twoI][twoJ]     +=  (bc[0][0] * bb1 + bc[0][2] * bb2) ;
                stiffness_local[twoI][twoJp1]   +=  (bc[0][1] * bb2 + bc[0][2] * bb1) ;
                stiffness_local[twoIp1][twoJ]   +=  (bc[1][0] * bb1 + bc[1][2] * bb2) ;
                stiffness_local[twoIp1][twoJp1] +=  (bc[1][1] * bb2 + bc[1][2] * bb1) ;
             }
          }

          //==============================================
          // CALCULATE  Kup, Ktt, Ktp
          //==============================================

          if(calcExtraMatrices)
          {
             for(jj=0;jj<sizeM;jj++)
             {
                fact1 = Nbar[jj] * dvol0;

                for(ii=0;ii<nlbf;ii++)
                {
                   twoI   = 2*ii;
                   twoIp1 = twoI+1;

                   Kup(twoI, jj)   += ( fact1 * dN_dx[ii] );
                   Kup(twoIp1, jj) += ( fact1 * dN_dy[ii] );
                }

                fact2 = BULK * fact1;

                for(ii=0;ii<sizeM;ii++)
                {
                   Ktt(ii,jj) += fact2 * Nbar[ii];
                   Ktp(ii,jj) -= fact1 * Nbar[ii];
                }
             }
          }

          //==============================================
          // CALCULATE residuals Rt and Rp
          //==============================================

          fact1 = (BULK*volstr1 - pres) * dvol0;
          fact2 = (volstr - volstr1) * dvol0;

          for(ii=0;ii<sizeM;ii++)
          {
             resi2[ii]       -= Nbar[ii] * fact1;
             resi2[sizeM+ii] -= Nbar[ii] * fact2;
          }

          count++;
          count1++;
          ll += nivGP;
      }
   }


  // fill the lower triangular entries(use of symmetry)
  for(ii=0;ii<nsize;ii++)
  {
     for(jj=ii+1;jj<nsize;jj++)
     {
        stiffness_local[jj][ii] = stiffness_local[ii][jj];
     }
  }

  calcExtraMatrices = false;

//printStiffnessMatrix();

//  printForceVector();

//  cout << endl;  cout << Kup << endl;  cout << endl;
*/
  return 0;
}


int NurbsElem2DStructMixed3field::calcStiffnessAndResidual2()
{
/*
  double BULK = matDat[0];

   int  err, isw, count, count1, ll, ii, jj, twoI, twoJ, twoIp1, twoJp1, index, gp1, gp2, mm;

   int  sizeM = surf2->nlbf;

   double  F[4], detF, F33, dvol0, dvol, dt, Jac, dummy, pres, bb1, bb2, bb3, detFbar, phat, pbar, r1d3=1.0/3.0, r2d3=2.0/3.0;

   double  strbar[4], strhat[4], strdev[4], bc[2][3], N[nlbf], dN_dx[nlbf], dN_dy[nlbf], Nbar[sizeM];

   double fact, fact1, fact2, fact3, fact4, utemp, vtemp;

   double  Idev[4][4], D12[4], D22, D11[4][4], cctmp[4][4], alpha, alpha1;

   Idev2D(Idev);

   if(calcExtraMatrices)
   {
      Kut.resize(nsize, sizeM);
      Kup.resize(nsize, sizeM);
      Ktt.resize(sizeM, sizeM);
      Ktp.resize(sizeM, sizeM);
   }

   for(ii=0;ii<nsize;ii++)
     stiffness_local[ii].zero();

   resi.zero();
   resi2.zero();
   Kut.setZero();
   Kup.setZero();
   Ktt.setZero();
   Ktp.setZero();


   count = 1;   ll = 0;   err = 0;   isw = 3;

   double *gaussweights = &(surf0->gaussweights[0]);

   count1 = 0;
   for(gp2=0;gp2<nGP2;gp2++)
   {
      for(gp1=0;gp1<nGP1;gp1++)
      {
          index = count1*2;

          utemp = knotsAtGPs[index];
          vtemp = knotsAtGPs[index+1];

          surf1->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, Jac);

          dvol = Jac * gaussweights[count1] * JacMultFact;

          surf0->deformationGradient(startindex[0], startindex[1], 0, dN_dx, dN_dy, F, detF);

          dvol0 = dvol/detF;

          surf2->ShapeFunctions(utemp, vtemp, Nbar);

          detFbar  = 1.0 + surf2->computeValue(1,utemp,vtemp);
          pres     = surf2->computeValue(2,utemp,vtemp);

          // calculate Fbar

          alpha = detFbar/detF;
          phat = pres/alpha;

          alpha1 = pow(alpha,r1d3);

//        printf(" alphas ");        printf("\t%12.8f\t%12.8f\n\n", alpha, alpha1);
//        printf(" dvol ");        printf("\t%12.8f\t%12.8f\n\n", dvol0, dvol);
          F33 = alpha1;
          for(ii=0;ii<4;ii++)
             F[ii] = alpha1*F[ii];

          dt = mpapTime.dt;
          matlib2d_(matDat, F, &F33, strbar, D11[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);
          if(err !=0)           return 1;

          dvol *= F33;

          pbar = (strbar[0]+strbar[1]+strbar[2])/3.0;

          strdev[0] = strbar[0] - pbar;
          strdev[1] = strbar[1] - pbar;
          strdev[2] = strbar[2] - pbar;
          strdev[3] = strbar[3];

          strhat[0] = strdev[0] + phat;
          strhat[1] = strdev[1] + phat;
          strhat[2] = strdev[2] + phat;
          strhat[3] = strbar[3];

//        printf(" stresses ");        printf("\t%12.8f\t%12.8f\t%12.8f\t%12.8f\n\n", stre[0], stre[1], stre[2], pres);
//        printf(" values   ");        printf("\t%12.8f\t%12.8f\t%12.8f\t%12.8f\t%12.8f\n\n", detF, detFbar, pbar, phat, dvol);

          // D22
          //--------------------------

          D22 = 0.0;
          for(ii=0;ii<3;ii++)
          {
             for(jj=0;jj<3;jj++)
                D22 += D11[ii][jj];
          }
          D22 /= 9.0;

//        printf(" D22 ");        printf("\t%12.8f\t%12.8f\t%12.8f\n\n", D22, pbar, BULK);

          D22 -= r1d3 * pbar;

//        printf(" D22 ");        printf("\t%12.8f\t%12.8f\t%12.8f\n\n", D22, pbar, BULK);

          // D11
          //--------------------------

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

          fact = 2.0 * (pbar - phat);
          fact1 = (r2d3*pbar - phat);

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

          // D12
          //--------------------------

          for(ii=0;ii<4;ii++)
          {
             D12[ii] = 0.0;
             for(jj=0;jj<3;jj++)
                D12[ii] += cctmp[ii][jj];
             D12[ii] *= r1d3;
          }
        //printf(" D12 ");        printf("\t%12.8f\t%12.8f\t%12.8f\t%12.8f\n\n", D12[0], D12[1], D12[2], D12[3]);
          for(ii=0;ii<4;ii++)
             D12[ii] += r2d3 * strdev[ii];
        //printf(" D12 ");        printf("\t%12.8f\t%12.8f\t%12.8f\t%12.8f\n\n", D12[0], D12[1], D12[2], D12[3]);


          // multiply with volume and corresponding factors of detF and detFbar

          fact = dvol0*detFbar;

          for(ii=0;ii<4;ii++)
          {
             strbar[ii] *= dvol;
             strhat[ii] *= fact;

             for(jj=0;jj<4;jj++)
               D11[ii][jj] *= fact;

             D12[ii] *= dvol0;
          }
          D22 *= (dvol0/detFbar);

        printf(" D22 ");        printf("\t%12.8f\t%12.8f\t%12.8f\t%12.8f\n\n", D22, pbar, BULK, detFbar);

          //==============================================
          // CALCULATE TANGENT STIFFNESS and RESIDUAL
          //==============================================

          for(ii=0;ii<nlbf;ii++)
          {
             bb1 = dN_dx[ii];
             bb2 = dN_dy[ii];

             bc[0][0] = (bb1 * D11[0][0] + bb2 * D11[3][0]);
             bc[0][1] = (bb1 * D11[0][1] + bb2 * D11[3][1]);
             bc[0][2] = (bb1 * D11[0][3] + bb2 * D11[3][3]);

             bc[1][0] = (bb2 * D11[1][0] + bb1 * D11[3][0]);
             bc[1][1] = (bb2 * D11[1][1] + bb1 * D11[3][1]);
             bc[1][2] = (bb2 * D11[1][3] + bb1 * D11[3][3]);

             twoI   = 2*ii;
             twoIp1 = twoI+1;

             fact1 = (bb1 * strhat[0] + bb2 * strhat[3]);
             fact2 = (bb1 * strhat[3] + bb2 * strhat[1]) ;

             resi[twoI]   -= (bb1 * strhat[0] + bb2 * strhat[3]);
             resi[twoIp1] -= (bb1 * strhat[3] + bb2 * strhat[1]) ;

             for(jj=0;jj<nlbf;jj++)
             {
                twoJ   = 2*jj;
                twoJp1 = twoJ+1;

                bb1 = dN_dx[jj];
                bb2 = dN_dy[jj];

                fact = fact1 * bb1 + fact2 * bb2;

                stiffness_local[twoI][twoJ]     +=  (bc[0][0] * bb1 + bc[0][2] * bb2 + fact) ;
                stiffness_local[twoI][twoJp1]   +=  (bc[0][1] * bb2 + bc[0][2] * bb1) ;
                stiffness_local[twoIp1][twoJ]   +=  (bc[1][0] * bb1 + bc[1][2] * bb2) ;
                stiffness_local[twoIp1][twoJp1] +=  (bc[1][1] * bb2 + bc[1][2] * bb1 + fact) ;
             }
          }

          //==============================================
          // CALCULATE Kut, Kup
          //==============================================

          fact = dvol0*detF;

          for(ii=0;ii<nlbf;ii++)
          {
             bb1 = dN_dx[ii];
             bb2 = dN_dy[ii];

             twoI   = 2*ii;
             twoIp1 = twoI+1;

             fact1 = (bb1 * D12[0] + bb2 * D12[3]) ;
             fact2 = (bb1 * D12[3] + bb2 * D12[1]) ;

             fact3 = bb1 * fact;
             fact4 = bb2 * fact;

             for(jj=0;jj<sizeM;jj++)
             {
                bb3 = Nbar[jj];

                Kut(twoI,jj)    +=  fact1 * bb3;
                Kut(twoIp1,jj)  +=  fact2 * bb3;

                Kup(twoI,jj)    +=  fact3 * bb3;
                Kup(twoIp1,jj)  +=  fact4 * bb3;
             }
          }

          //==============================================
          // CALCULATE Ktt, Ktp
          //==============================================

          for(ii=0;ii<sizeM;ii++)
          {
             bb1 = Nbar[ii];
             fact1 = bb1 * D22;
             fact2 = bb1 * dvol0;

             for(jj=0;jj<sizeM;jj++)
             {
                bb2 = Nbar[jj];
                Ktt(ii,jj) += fact1 * bb2;
                Ktp(ii,jj) -= fact2 * bb2;
             }
          }

          //==============================================
          // CALCULATE residuals Rt and Rp
          //==============================================

          fact1 = (pbar - pres) * dvol0;
          fact2 = (detF - detFbar) * dvol0;

          for(ii=0;ii<sizeM;ii++)
          {
             bb1 = Nbar[ii];
             resi2[ii]       -= bb1 * fact1;
             resi2[sizeM+ii] -= bb1 * fact2;
          }

          count++;
          count1++;
          ll += nivGP;
      }
   }

  calcExtraMatrices = false;

printStiffnessMatrix();

//  printForceVector();

  cout << " Kup " << endl;
  cout << endl;  cout << Kup << endl;  cout << endl;
  cout << " Kut " << endl;
  cout << endl;  cout << Kut << endl;  cout << endl;
  cout << " Ktt " << endl;
  cout << endl;  cout << Ktt << endl;  cout << endl;
  cout << " Ktp " << endl;
  cout << endl;  cout << Ktp << endl;  cout << endl;

//cout << resi2 << endl;
*/
  return 0;
}



int NurbsElem2DStructMixed3field::calcInternalForces()
{

  return 0;
}




int NurbsElem2DStructMixed3field::calcOutput(double u1, double v1)
{

  return 0;
}







void NurbsElem2DStructMixed3field::discreteContourplot(int vartype, int varindex, int index, int nCol, double umin, double umax)
{
  if(index > nivGP)
  {
     cout << '\t' << " Error in NurbsElem2DStructMixed3field::contourplot " << endl;
     return;
  }

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








void NurbsElem2DStructMixed3field::projectToKnots(bool extrapolateFlag, int vartype, int varindex, int index)
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

//   cout << '\t' << vals2project << endl; cout << endl;

  return;
}




void NurbsElem2DStructMixed3field::projectStress(int varindex, double* outval)
{
/*
   int nivEL = nGP * nivGP;
   for(int ii=0;ii<nivEL;ii++)
     intVar2[ii] = intVar1[ii];

   int  err, isw, count, count1, ll, index, gp1, gp2, ii;

   int  sizeM = surf2->nlbf;

   double  F[4], detF, F33, dvol0, dvol, dt, Jac, pres, detFbar, phat, alpha, alpha1, r1d3=1.0/3.0;

   double  cc[4][4], stre[4], N[nlbf], dN_dx[nlbf], dN_dy[nlbf], Nbar[sizeM];

   double fact, utemp, vtemp, volstr, volstr1, dummy, pbar;


   count = 1;   ll = 0;   err = 0;   isw = 3;

   count1 = 0;

   if(finite)
   {
      // loop over Gauss points
      for(gp2=0;gp2<nGP2;gp2++)
      {
         for(gp1=0;gp1<nGP1;gp1++)
         {
             index = count1*2;

             utemp = knotsAtGPs[index];
             vtemp = knotsAtGPs[index+1];

             surf0->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, Jac);

             surf1->deformationGradient(startindex[0], startindex[1], 1, dN_dx, dN_dy, F, detF);

             detFbar  = 1.0 + surf2->computeValue(1,utemp,vtemp);
             pres     = surf2->computeValue(2,utemp,vtemp);

             // calculate Fbar

             alpha = detFbar/detF;
             phat = pres/alpha;

             alpha1 = pow(alpha,r1d3);

             F33 = alpha1;
             for(ii=0;ii<4;ii++)
               F[ii] = alpha1*F[ii];

             dt = mpapTime.dt;
             matlib2d_(matDat, F, &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);

             pbar = (stre[0]+stre[1]+stre[2])/3.0;

             fact = phat - pbar;

             stre[0] += fact;
             stre[1] += fact;
             stre[2] += fact;

             if(varindex < 4)
                outval[count1] = stre[varindex];
             else if(varindex == 4)
                outval[count1] = sqrt((pow(stre[0]-stre[1],2.0) + pow(stre[1]-stre[2], 2.0) + pow(stre[2]-stre[0], 2.0) + 6* stre[3]*stre[3])/2);
             else if(varindex == 5)
                outval[count1] = pres;
             else
             {
                cout << '\t' << "    NurbsElem2DStructMixed3field::projectStress .... : Error in 'varindex' " << endl;
               return;
             }

             count++;
             count1++;
             ll += nivGP;
         }
      }
   }
   else
   {
      // loop over Gauss points
      for(gp2=0;gp2<nGP2;gp2++)
      {
         for(gp1=0;gp1<nGP1;gp1++)
         {
             index = count1*2;

             utemp = knotsAtGPs[index];
             vtemp = knotsAtGPs[index+1];

             surf0->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, Jac);

             surf1->deformationGradient(startindex[0], startindex[1], 1, dN_dx, dN_dy, F, detF);

             volstr1  = surf2->computeValue(1,utemp,vtemp);
             pres     = surf2->computeValue(2,utemp,vtemp);

             // modify F

             volstr = (F[0]+F[3]-2.0);

             fact = (volstr1 - volstr)/3.0;

             F[0] += fact;
             F[3] += fact;
             F33 = 1.0 + fact;

             dt = mpapTime.dt;
             matlib2d_(matDat, F, &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, (Element*) this);

             pbar = (stre[0]+stre[1]+stre[2])/3.0;

             fact = pres - pbar;

             stre[0] += fact;
             stre[1] += fact;
             stre[2] += fact;

             if(varindex < 4)
                outval[count1] = stre[varindex];
             else if(varindex == 4)
                outval[count1] = sqrt((pow(stre[0]-stre[1],2.0) + pow(stre[1]-stre[2], 2.0) + pow(stre[2]-stre[0], 2.0) + 6* stre[3]*stre[3])/2);
             else if(varindex == 5)
                outval[count1] = pres;
             else
             {
                cout << '\t' << "    NurbsElem2DStructMixed3field::projectStress .... : Error in 'varindex' " << endl;
               return;
             }

             count++;
             count1++;
             ll += nivGP;
         }
      }
   }
*/
  return;
}




void NurbsElem2DStructMixed3field::projectStrain(int vartype, int varindex, double* outval)
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










void NurbsElem2DStructMixed3field::projectIntVar(int index, double* outval)
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




void NurbsElem2DStructMixed3field::AssembleElementMatrix2(int index, MatrixXd& stiff, MatrixXd& CC)
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



void NurbsElem2DStructMixed3field::AssembleElementVector2(bool firstIter, int flag, VectorXd& rhs, double* reac, VectorXd& uu)
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
      cerr << " ERROR in flag in 'NurbsElem2DStructMixed3field::AssembleElementVector2' " << endl;
      return;
   }


//  cout << " resi " << resi << endl;

  return;
}






void NurbsElem2DStructMixed3field::AssembleElementMatrix(int index, MatrixSparseArray<double>& mtx)
{
     int nn=0, aa, bb, size2;

     size2 = surf2->nsize;

     for(aa=0;aa<nsize;aa++)
     {
         for(bb=0;bb<nsize;bb++)
         {
             nn = forassembly[aa][bb];
             if(nn != -1)
                mtx.x[nn-1] += stiffness_local[aa][bb];
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
                mtx.x[nn-1] += Ktt(aa,bb);

             nn = forassyKtp[aa][bb];
             if(nn != -1)
                mtx.x[nn-1] += Ktp(aa,bb);

             nn = forassyKpt[bb][aa];
             if(nn != -1)
                mtx.x[nn-1] += Ktp(aa,bb);
         }
     }

     if(finite)
     {
         for(aa=0;aa<nsize;aa++)
         {
            for(bb=0;bb<size2;bb++)
            {
                nn = forassyKut[aa][bb];
                if(nn != -1)
                   mtx.x[nn-1] += Kut(aa,bb);

                nn = forassyKtu[bb][aa];
                if(nn != -1)
                   mtx.x[nn-1] += Kut(aa,bb);
            }
         }
     }

  return;
}





void NurbsElem2DStructMixed3field::AssembleElementVector(bool firstIter, bool flag, double* rhs, double* reac, int start1, int start2)
{
   // flag == true  ---> just external force vector
   // flag == false ---> internal load vector + contributions from nodes with specified displacement BCs

   int *tt;

   tt = &(surf0->LM[elenum][0]);

   if(flag)
   {
      for(int aa=0;aa<nsize;aa++)
      {
         if(tt[aa] != -1)
            rhs[tt[aa]] += resi[aa];
      }
   }
   else
   {
      double fact;
      int aa, bb, *tt2, size2, ind;

      tt2 = &(surf2->LM[elenum][0]);
      size2 = surf2->nsize;

      for(aa=0;aa<nsize;aa++)
      {
         if(tt[aa] != -1)
           rhs[tt[aa]] += resi[aa];

         // add up reaction forces
         reac[forassy[aa]] += resi[aa];
      }

      for(aa=0;aa<size2;aa++)
      {
         if(tt2[aa] != -1)
         {
            rhs[tt2[aa]+start1] += resi2[aa];
            rhs[tt2[aa]+start2] += resi2[aa+size2];
         }
      }

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
                      rhs[tt[bb]] -= stiffness_local[bb][aa] * fact;
                 }
             }
         }

         for(aa=0;aa<size2;aa++)
         {
             ind = tt2[aa] + start2;

             for(bb=0;bb<nsize;bb++)
             {
                 if(tt[bb] == -1)
                   rhs[ind] -= Kup(bb,aa) * mpapTime.dt * primvar[bb];
             }
         }
      }
      tt2 = NULL;
   }

//  cout << " resi " << resi << endl;

  tt= NULL;

  return;
}




void NurbsElem2DStructMixed3field::toPostprocess(int vartype, int varindex, int type, SparseMatrixXd&  coeffMat, VectorXd& rhsVec)
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




