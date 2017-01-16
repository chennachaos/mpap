
#include <math.h>
#include "Debug.h"
#include "Plot.h"
#include "FunctionsElement.h"
#include "MpapTime.h"
#include "NurbsElem2DStructMixed2field.h"
#include "NurbsShapeFunctions.h"
#include <assert.h>
#include "ComputerTime.h"

#include "util.h"

#include "ExactSolutionsElasticity.h"

using namespace std;

extern Plot     plot;
extern MpapTime mpapTime;
extern ComputerTime       computerTime;



NurbsElem2DStructMixed2field::NurbsElem2DStructMixed2field()
{
  if (debug) cout << " constructor NurbsElem2DStructMixed2field\n\n";

//   cout << " constructor NurbsElem2DStructMixed2field\n\n";

  calcExtraMatrices = true;
}



NurbsElem2DStructMixed2field::~NurbsElem2DStructMixed2field()
{
  if (debug) cout << " destructor NurbsElem2DStructMixed2field\n\n";

  // cout << " destructor NurbsElem2DStructMixed2field\n\n";

}

int NurbsElem2DStructMixed2field::calcStiffnessMatrix(double dt)
{

  return 0;
}



void NurbsElem2DStructMixed2field::contourplot(int index, int nCol, double umin, double umax)
{
  return;
}

int NurbsElem2DStructMixed2field::calcMassMatrix(int lumpInd, double dt)
{


  return 0;
}


int NurbsElem2DStructMixed2field::calcStiffnessAndResidual()
{
  if(axsy)
  {
    //cout << " AXISYMMETRIC case is not supported by NurbsElem2DStructMixed2field  " << endl;

    if(finite)
      NurbsElem2DStructMixed2field::calcStiffnessAndResidualAxsyFS();
    else
      NurbsElem2DStructMixed2field::calcStiffnessAndResidualAxsySS();
  }
  else
  {
    if(finite)
      NurbsElem2DStructMixed2field::calcStiffnessAndResidual2();
    else
      NurbsElem2DStructMixed2field::calcStiffnessAndResidual1();
  }

  return 0;
}





//
int NurbsElem2DStructMixed2field::calcStiffnessAndResidual1()
{
  // elasticity mixed formulation

  int  err, isw, count, count1, ll, ii, jj, twoI, twoJ, twoIp1, twoJp1, index, gp1, gp2, mm;
  int  sizep = surf2->nlbf;

  double  F[4], detF, F33, fact, dvol0, dt, Jac, dummy, pres, bb1, bb2, bb3, volstr, force[2];
  double  cc[4][4], stre[4], bc[2][3], Idev[4][4], cctmp[4][4];
  double  eps, BULK = matDat[0];

  vector<double>  N(nlbf), dN_dx(nlbf), dN_dy(nlbf), Nbar(sizep);

  Idev2D(Idev);

  double *gaussweights = &(surf0->gaussweights[0]);

  //double  mu,  E=210.0, nu=0.49999;

  //BULK  = E/3.0/(1.0-2.0*nu);
  //mu    = E/2.0/(1.0+nu);

  //matDat[0] = BULK;
  //matDat[1] = mu;
  //matDat[2] = nu;
  
  //ElasticityEx1  analy(matDat[1]);
  //ElasticityEx2  analy(matDat[0], matDat[1]);
  
  EPOINT  EP;
  
  //BULK = matDat[0];
  //if(elmDat[9] == 1)
    eps = 1.0/BULK;
  //else
    //eps = 0.0;

   Klocal.setZero();
   Flocal.setZero();
   resi2.zero();

   Kup.resize(nsize, sizep);
   Kup.setZero();

   Kpp.resize(sizep, sizep);
   Kpp.setZero();

   count = 1;   ll = 0;   err = 0;   isw = 3;
   dt = mpapTime.dt;

  //cout << nGP1 << '\t' << nGP2 << endl;
  //cout << rho0 << '\t' << bforce[0] << '\t' << bforce[1] << endl;

   count1 = 0;
   for(gp2=0;gp2<nGP2;gp2++)
   {
      for(gp1=0;gp1<nGP1;gp1++)
      {
          index = count1*2;

          //cout << gp1 << '\t' << gp2 << endl;

          surf0->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), &N[0], &dN_dx[0], &dN_dy[0], Jac);

          dvol0 = Jac * gaussweights[count1] * JacMultFact;

          EP = surf0->SurfacePoint(knotsAtGPs[index], knotsAtGPs[index+1]).CalcEuclid();

        //force[0] = analy.forceX(EP.x, EP.y);
        //force[1] = analy.forceY(EP.x, EP.y);

          //cout << " lllllllll " << endl;

          surf1->deformationGradient(startindex[0], startindex[1], 1, &dN_dx[0], &dN_dy[0], F, detF);

          //cout << " lllllllll " << endl;

          volstr = (F[0]+F[3]-2.0);

          F33 = 1.0;

          matlib2d_(matDat, F, &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);
          if(err !=0)           return 1;

          pres = surf2->computeValueAndShanpeFns(1, knotsAtGPs[index], knotsAtGPs[index+1], &Nbar[0]);

          dummy = pres - (stre[0]+stre[1]+stre[2])/3.0 ;

          stre[0] += dummy;
          stre[1] += dummy;
          stre[2] += dummy;

          //printf(" stresses ");        printf("\t%12.8f\t%12.8f\t%12.8f\t%12.8f\n\n", stre[0], stre[1], stre[2], pres);

          //printf("\t%12.8f\t%12.8f\t%12.8f\t%12.8f\n\n", Jac, JacMultFact, dvol0, volstr);

          //for(ii=0;ii<nlbf;ii++)
            //printf("\t%12.8f\t%12.8f\t%12.8f\n\n", N[ii], dN_dx[ii], dN_dy[ii]);

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

          for(ii=0;ii<4;ii++)
          {
             stre[ii] *= dvol0;
             for(jj=0;jj<4;jj++)
               cc[ii][jj] *= dvol0;
          }

          //==============================================
          // CALCULATE TANGENT STIFFNESS and RESIDUAL
          //==============================================

          fact = rho0 * dvol0;
          //force[0] = bforce[0] * fact ;
          //force[1] = bforce[1] * fact ;
          //force[0] = 0.0;
          //force[1] = 0.0;

          for(ii=0;ii<nlbf;ii++)
          {
             bb1 = dN_dx[ii];
             bb2 = dN_dy[ii];
             bb3 = N[ii]*dvol0;

             bc[0][0] = (bb1 * cc[0][0] + bb2 * cc[3][0]);
             bc[0][1] = (bb1 * cc[0][1] + bb2 * cc[3][1]);
             bc[0][2] = (bb1 * cc[0][3] + bb2 * cc[3][3]);

             bc[1][0] = (bb2 * cc[1][0] + bb1 * cc[3][0]);
             bc[1][1] = (bb2 * cc[1][1] + bb1 * cc[3][1]);
             bc[1][2] = (bb2 * cc[1][3] + bb1 * cc[3][3]);

             twoI   = 2*ii;
             twoIp1 = twoI+1;

             Flocal[twoI]   += (bb3*force[0] - bb1*stre[0] - bb2*stre[3]) ;
             Flocal[twoIp1] += (bb3*force[1] - bb1*stre[3] - bb2*stre[1]) ;

             for(jj=0;jj<nlbf;jj++)
             {
                twoJ   = 2*jj;
                twoJp1 = twoJ+1;

                bb1 = dN_dx[jj];
                bb2 = dN_dy[jj];

                Klocal(twoI,twoJ)     +=  (bc[0][0] * bb1 + bc[0][2] * bb2) ;
                Klocal(twoI,twoJp1)   +=  (bc[0][1] * bb2 + bc[0][2] * bb1) ;
                Klocal(twoIp1,twoJ)   +=  (bc[1][0] * bb1 + bc[1][2] * bb2) ;
                Klocal(twoIp1,twoJp1) +=  (bc[1][1] * bb2 + bc[1][2] * bb1) ;
             }
             // compute Kup and Kpp matrices

             for(jj=0;jj<sizep;jj++)
             {
                fact = Nbar[jj] * dvol0;

                Kup(twoI, jj)   += ( fact * dN_dx[ii] );
                Kup(twoIp1, jj) += ( fact * dN_dy[ii] );
             }
          }

          bb1 = dvol0*eps;
          bb2 = (volstr - pres*eps)*dvol0;

          for(ii=0;ii<sizep;ii++)
          {
            resi2[ii] -= (bb2 * Nbar[ii]);

            fact = Nbar[ii] * bb1;

            for(jj=0;jj<sizep;jj++)
              Kpp(ii,jj) -= ( fact * Nbar[jj] );
          }

          count++;
          count1++;
          ll += nivGP;
     }
   }

  Kpu = Kup.transpose();

  //  printMatrix(Klocal);  printf("\n\n");  printVector(Flocal);

  return 0;
}
//




int NurbsElem2DStructMixed2field::toComputeInfSupCondition()
{
   // to compute the inf-sup constant
   
   int  err, isw, count, count1, ll, ii, jj, twoI, twoJ, twoIp1, twoJp1, index, gp1, gp2, mm;
   int  sizep = surf2->nlbf;

   double  F[4], detF, F33, fact, dvol0, dt, Jac, dummy, pres, b1, b2, b3, b4, b5, volstr;
   double  eta = elmDat[7];

   vector<double>  N(nlbf), dN_dx(nlbf), dN_dy(nlbf), Nbar(sizep);

   double *gaussweights = &(surf0->gaussweights[0]);

   Kup.resize(nsize, sizep);
   Kpp.resize(sizep, sizep);

   Klocal.setZero();
   Kup.setZero();
   Kpp.setZero();
   Flocal.setZero();
   resi2.zero();

   count = 1;   ll = 0;   err = 0;   isw = 3;
   dt = mpapTime.dt;

   count1 = 0;
   for(gp2=0;gp2<nGP2;gp2++)
   {
      for(gp1=0;gp1<nGP1;gp1++)
      {
          index = count1*2;

          surf0->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), &N[0], &dN_dx[0], &dN_dy[0], Jac);
          dvol0 = Jac * gaussweights[count1] * JacMultFact;

          //surf0->ShapeFunctions(knotsAtGPs[index], knotsAtGPs[index+1], N);
          surf2->ShapeFunctions(knotsAtGPs[index], knotsAtGPs[index+1], &Nbar[0]);

          //==============================================
          // CALCULATE TANGENT STIFFNESS and RESIDUAL
          //==============================================

          for(ii=0;ii<nlbf;ii++)
          {
             twoI   = 2*ii;
             twoIp1 = twoI+1;

             b1 = dN_dx[ii]*dvol0;
             b2 = dN_dy[ii]*dvol0;
             b3 = eta*N[ii]*dvol0;

             for(jj=0;jj<nlbf;jj++)
             {
                twoJ   = 2*jj;
                twoJp1 = twoJ+1;

                fact = b1*dN_dx[jj] + b2*dN_dy[jj] + b3*N[jj];

                Klocal(twoI,twoJ)     += fact;
                Klocal(twoI,twoJp1)   += 0.0;
                Klocal(twoIp1,twoJ)   += 0.0;
                Klocal(twoIp1,twoJp1) += fact;
             }
             // compute Kup and Kpp matrices

             for(jj=0;jj<sizep;jj++)
             {
                Kup(twoI, jj)   += ( b1 * Nbar[jj] );
                Kup(twoIp1, jj) += ( b2 * Nbar[jj] );
             }
          }
             
          for(ii=0;ii<sizep;ii++)
          {
             fact = Nbar[ii] * dvol0;

             for(jj=0;jj<sizep;jj++)
                Kpp(ii,jj) += ( fact * Nbar[jj] );
          }

          count1++;
     }
   }

//  printStiffnessMatrix();
//  printf("\n\n");
//  printForceVector();

  return 0;
}





//
int NurbsElem2DStructMixed2field::calcStiffnessAndResidual2()
{
  // mixed formulation - no stabilization
  
  double BULK = matDat[0];

   int  err, isw, count, count1, ll, ii, jj, twoI, twoJ, twoIp1, twoJp1, index, gp1, gp2, mm;

   int  sizep = surf2->nlbf;

   double  F[4], detF, F33, dvol0, dvol, dt, Jac, dummy, pres, bb1, bb2, bb3, pbar, r2d3 = 2.0/3.0;
   double  stre[4], strdev[4], bc[2][3];
   double  fact, fact1, fact2, fact3, fact4, eps;
   double  Idev[4][4], D11[4][4], cctmp[4][4];

   vector<double>  N(nlbf), dN_dx(nlbf), dN_dy(nlbf), Nbar(sizep);

   //if(elmDat[9] == 1)
     eps = 1.0/BULK;
   //else
     //eps = 0.0;

   Idev2D(Idev);

   double *gaussweights = &(surf0->gaussweights[0]);

   Kup.resize(nsize, sizep);
   Kpp.resize(sizep, sizep);

   Klocal.setZero();
   Kup.setZero();
   Kpp.setZero();
   Flocal.setZero();
   resi2.zero();

   count = 1;   ll = 0;   err = 0;   isw = 3;
   dt = mpapTime.dt;

   count1 = 0;
   for(gp2=0;gp2<nGP2;gp2++)
   {
      for(gp1=0;gp1<nGP1;gp1++)
      {
          index = count1*2;

          //surf1->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), dN_dx, dN_dy, Jac);

          //dvol = Jac * gaussweights[count1] * JacMultFact;

          //surf0->deformationGradient(startindex[0], startindex[1], 0, dN_dx, dN_dy, F, detF);

          //dvol0 = dvol/detF;

          //cout << " llllllllll " << endl;
          surf0->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), &N[0], &dN_dx[0], &dN_dy[0], Jac);

          dvol0 = Jac * gaussweights[count1] * thick * JacMultFact;

          surf1->deformationGradient(startindex[0], startindex[1], 1, &dN_dx[0], &dN_dy[0], F, detF);

          surf1->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), &N[0], &dN_dx[0], &dN_dy[0], Jac);

          dvol = Jac * gaussweights[count1] * thick * JacMultFact;

          pres = surf2->computeValueAndShanpeFns(1, knotsAtGPs[index], knotsAtGPs[index+1], &Nbar[0]);
          //cout << " llllllllll " << endl;

          //surf2->ShapeFunctions(utemp, vtemp, Nbar);
          //pres  = surf2->computeValue(1, utemp, vtemp);

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

          for(ii=0;ii<4;ii++)
          {
             stre[ii] *= dvol;
             for(jj=0;jj<4;jj++)
               D11[ii][jj] *= dvol;
          }

          //==============================================
          // CALCULATE STIFFNESS matrices and RESIDUALs
          //==============================================

          for(ii=0;ii<nlbf;ii++)
          {
             bb1 = dN_dx[ii];
             bb2 = dN_dy[ii];

             twoI   = 2*ii;
             twoIp1 = twoI+1;

             bc[0][0] = (bb1 * D11[0][0] + bb2 * D11[3][0]);
             bc[0][1] = (bb1 * D11[0][1] + bb2 * D11[3][1]);
             bc[0][2] = (bb1 * D11[0][3] + bb2 * D11[3][3]);

             bc[1][0] = (bb2 * D11[1][0] + bb1 * D11[3][0]);
             bc[1][1] = (bb2 * D11[1][1] + bb1 * D11[3][1]);
             bc[1][2] = (bb2 * D11[1][3] + bb1 * D11[3][3]);

             fact1 = (bb1 * stre[0] + bb2 * stre[3]) ;
             fact2 = (bb1 * stre[3] + bb2 * stre[1]) ;

             fact3 = bb1 * dvol;
             fact4 = bb2 * dvol;

             Flocal[twoI]   -= fact1 ;
             Flocal[twoIp1] -= fact2 ;

             for(jj=0;jj<nlbf;jj++)
             {
                twoJ   = 2*jj;
                twoJp1 = twoJ+1;

                bb1 = dN_dx[jj];
                bb2 = dN_dy[jj];

                fact = fact1 * bb1 + fact2 * bb2;

                Klocal(twoI,twoJ)     +=  (bc[0][0] * bb1 + bc[0][2] * bb2+ fact) ;
                Klocal(twoI,twoJp1)   +=  (bc[0][1] * bb2 + bc[0][2] * bb1) ;
                Klocal(twoIp1,twoJ)   +=  (bc[1][0] * bb1 + bc[1][2] * bb2) ;
                Klocal(twoIp1,twoJp1) +=  (bc[1][1] * bb2 + bc[1][2] * bb1+ fact) ;
             }

             for(jj=0;jj<sizep;jj++)
             {
                bb3 = Nbar[jj];

                Kup(twoI,jj)    +=  fact3 * bb3;
                Kup(twoIp1,jj)  +=  fact4 * bb3;
             }
          }

          fact = dvol0*eps;

          fact1 = (detF - 1.0 - pres*eps) * dvol0;

          //fact1 = (0.5*(detF - 1.0/detF) - pres/BULK) * dvol0;

          for(ii=0;ii<sizep;ii++)
          {
              fact2 = Nbar[ii] * fact;

              for(jj=0;jj<sizep;jj++)
                 Kpp(ii,jj) -= ( fact2 * Nbar[jj] );

              resi2[ii] -= Nbar[ii] * fact1;
          }

          count++;
          count1++;
          ll += nivGP;
      }
   }
   
   Kpu = Kup.transpose();

  return 0;
}
//




int NurbsElem2DStructMixed2field::calcStiffnessAndResidualAxsySS()
{

  return 0;
}



int NurbsElem2DStructMixed2field::calcStiffnessAndResidualAxsyFS()
{
   double BULK = matDat[0];

   int  err, isw, count, count1, ll, ii, jj, twoI, twoJ, twoIp1, twoJp1, index, gp1, gp2, mm;

   int  sizep = surf2->nlbf;

   double  F[4], detF, F33, dvol0, dvol, dt, Jac, dummy, pres, bb1, bb2, bb3, pbar, r2d3 = 2.0/3.0;
   double  stre[4], strdev[4], bc[2][4];
   double fact, fact1, fact2, fact3, fact4, fact5, fact6, fact7, rad0, rad, r1drad0, r1drad, utemp, vtemp;
   double  Idev[4][4], D11[4][4], cctmp[4][4];

  vector<double>  N(nlbf), dN_dx(nlbf), dN_dy(nlbf), Nbar(sizep);

   Idev2D(Idev);

   if(calcExtraMatrices)
   {
      Kup.resize(nsize, sizep);
      Kpp.resize(sizep, sizep);
   }

   double *gaussweights = &(surf0->gaussweights[0]);

   for(ii=0;ii<nsize;ii++)
     stiffness_local[ii].zero();

   Kup.setZero();
   Kpp.setZero();
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

          surf0->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), &N[0], &dN_dx[0], &dN_dy[0], Jac);

          fact = twoPI * gaussweights[count1] * JacMultFact;

          rad0 = surf0->SurfacePoint(utemp, vtemp).CalcEuclid().x;

          dvol0  = rad0 * fact * Jac;

          surf1->deformationGradient(startindex[0], startindex[1], 1, &dN_dx[0], &dN_dy[0], F, detF);

          surf1->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), &N[0], &dN_dx[0], &dN_dy[0], Jac);

          rad  = surf1->SurfacePoint(utemp, vtemp).CalcEuclid().x;

          dvol  = rad0 * fact * Jac;
          
          r1drad  = 1.0/rad;

          surf0->ShapeFunctions(utemp, vtemp, &N[0]);

          pres = surf2->computeValueAndShanpeFns(1, utemp, vtemp, &Nbar[0]);

          F33 = rad/rad0;

          matlib2d_(matDat, F, &F33, stre, D11[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);
          if(err !=0)           return 1;

          dvol *= F33;

          pbar = (stre[0]+stre[1]+stre[2])/3.0;

          strdev[0] = stre[0] - pbar;
          strdev[1] = stre[1] - pbar;
          strdev[2] = stre[2] - pbar;
          strdev[3] = stre[3];

          stre[0] = strdev[0] + pres;
          stre[1] = strdev[1] + pres;
          stre[2] = strdev[2] + pres;

//        printf(" stresses ");        printf("\t%12.8f\t%12.8f\t%12.8f\t%12.8f\n\n", stre[0], stre[1], stre[2], pres);
//        printf(" volumes   ");        printf("\t%12.8f\t%12.8f\t%12.8f\t%12.8f\n\n", F33, detF, dvol0, dvol);

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

          for(ii=0;ii<4;ii++)
          {
             stre[ii] *= dvol;
             for(jj=0;jj<4;jj++)
               D11[ii][jj] *= dvol;
          }

          //==============================================
          // CALCULATE STIFFNESS matrices and RESIDUALs
          //==============================================
          for(ii=0;ii<nlbf;ii++)
            N[ii] *= r1drad;

          for(ii=0;ii<nlbf;ii++)
          {
             bb1 = dN_dx[ii];
             bb2 = dN_dy[ii];
             fact = N[ii];

             bc[0][0] = (bb1 * D11[0][0] + bb2 * D11[3][0] + fact * D11[2][0]) ;
             bc[0][1] = (bb1 * D11[0][1] + bb2 * D11[3][1] + fact * D11[2][1]);
             bc[0][2] = (bb1 * D11[0][2] + bb2 * D11[3][2] + fact * D11[2][2]);
             bc[0][3] = (bb1 * D11[0][3] + bb2 * D11[3][3] + fact * D11[2][3]);

             bc[1][0] = (bb2 * D11[1][0] + bb1 * D11[3][0]);
             bc[1][1] = (bb2 * D11[1][1] + bb1 * D11[3][1]);
             bc[1][2] = (bb2 * D11[1][2] + bb1 * D11[3][2]);
             bc[1][3] = (bb2 * D11[1][3] + bb1 * D11[3][3]);

             twoI   = 2*ii;
             twoIp1 = twoI+1;

             fact1 = (bb1 * stre[0] + bb2 * stre[3]) ;
             fact2 = (bb1 * stre[3] + bb2 * stre[1]) ;

             fact3 = (bb1+fact) * dvol;
             fact4 = bb2 * dvol;

             fact6 = fact*stre[2];

             resi[twoI]   -= (fact1  + fact*stre[2]);
             resi[twoIp1] -= fact2 ;

             for(jj=0;jj<nlbf;jj++)
             {
                twoJ   = 2*jj;
                twoJp1 = twoJ+1;

                bb1 = dN_dx[jj];
                bb2 = dN_dy[jj];
                fact = N[jj];

                fact7 = fact1 * bb1 + fact2 * bb2;

                stiffness_local[twoI][twoJ]     +=  (bc[0][0] * bb1 + bc[0][3] * bb2 + bc[0][2] * fact + fact7 + fact6 * fact) ;
                stiffness_local[twoI][twoJp1]   +=  (bc[0][1] * bb2 + bc[0][3] * bb1) ;
                stiffness_local[twoIp1][twoJ]   +=  (bc[1][0] * bb1 + bc[1][3] * bb2 + bc[1][2] * fact) ;
                stiffness_local[twoIp1][twoJp1] +=  (bc[1][1] * bb2 + bc[1][3] * bb1 + fact7 ) ;
             }

             for(jj=0;jj<sizep;jj++)
             {
                bb3 = Nbar[jj];

                Kup(twoI,jj)    +=  fact3 * bb3;
                Kup(twoIp1,jj)  +=  fact4 * bb3;
             }
          }

          fact = dvol0/BULK;
          
          detF *= F33;

          fact1 = (detF - 1.0 - pres/BULK) * dvol0;

//          printf(" detF ");        printf("\t%12.8f\t%12.8f\t%12.8f\t%12.8f\n", detF,pres,dvol0,fact1);
          
          for(ii=0;ii<sizep;ii++)
          {
              fact2 = Nbar[ii] * fact;
              
              for(jj=0;jj<sizep;jj++)
                 Kpp(ii,jj) -= ( fact2 * Nbar[jj] );

              resi2[ii] -= Nbar[ii] * fact1;
          }

          count++;
          count1++;
          ll += nivGP;
      }
   }

   calcExtraMatrices = false;

//  printStiffnessMatrix();

//  printForceVector();

//  cout << endl;  cout << Kup << endl;  cout << endl;

//  cout << resi2 << endl;

  return 0;
}


int NurbsElem2DStructMixed2field::calcInternalForces()
{
   double BULK = matDat[0];

   int  err, isw, count, count1, ll, ii, jj, twoI, twoIp1, index, gp1, gp2;
   int  sizep = surf2->nlbf;

   double  F[4], detF, F33, dvol0, dvol, dt, Jac, dummy, pres, bb1, bb2, bb3, pbar, r2d3 = 2.0/3.0;
   double  stre[4], strdev[4];
   double  fact, fact1, fact2, D11[4][4], eps;

   vector<double>  N(nlbf), dN_dx(nlbf), dN_dy(nlbf), Nbar(sizep);

   double *gaussweights = &(surf0->gaussweights[0]);

   Flocal.setZero();
   resi2.zero();

   //if(elmDat[9] == 1)
     eps = 1.0/BULK;
   //else
     //eps = 0.0;

   count = 1;   ll = 0;   err = 0;   isw = 3;
   dt = mpapTime.dt;

   count1 = 0;
   for(gp2=0;gp2<nGP2;gp2++)
   {
      for(gp1=0;gp1<nGP1;gp1++)
      {
          index = count1*2;

          surf0->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), &N[0], &dN_dx[0], &dN_dy[0], Jac);

          dvol0 = Jac * gaussweights[count1] * thick * JacMultFact;

          surf1->deformationGradient(startindex[0], startindex[1], 1, &dN_dx[0], &dN_dy[0], F, detF);

          surf1->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), &N[0], &dN_dx[0], &dN_dy[0], Jac);

          dvol = Jac * gaussweights[count1] * thick * JacMultFact;

          pres = surf2->computeValueAndShanpeFns(1, knotsAtGPs[index], knotsAtGPs[index+1], &Nbar[0]);

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

//        printf(" stresses ");        printf("\t%12.8f\t%12.8f\t%12.8f\t%12.8f\n\n", stre[0], stre[1], stre[2], pres);

          for(ii=0;ii<4;ii++)
             stre[ii] *= dvol;

          for(ii=0;ii<nlbf;ii++)
          {
             bb1 = dN_dx[ii];
             bb2 = dN_dy[ii];

             twoI   = 2*ii;
             twoIp1 = twoI+1;

             Flocal[twoI]   -= (bb1 * stre[0] + bb2 * stre[3]);
             Flocal[twoIp1] -= (bb1 * stre[3] + bb2 * stre[1]) ;
          }

          fact1 = (detF - 1.0 - pres*eps) * dvol0;

          for(ii=0;ii<sizep;ii++)
             resi2[ii] -= Nbar[ii] * fact1;

          count++;
          count1++;
          ll += nivGP;
      }
   }
/*
  printForceVector();
  cout << endl;
  cout << endl;
  cout << resi2 << endl;
  cout << endl;
  cout << endl;
*/
  return 0;
}




int NurbsElem2DStructMixed2field::calcOutput(double u1, double v1)
{

  return 0;
}







void NurbsElem2DStructMixed2field::discreteContourplot(int vartype, int varindex, int index, int nCol, double umin, double umax)
{
  if(index > nivGP)
  {
     cout << '\t' << " Error in NurbsElem2DStructMixed2field::contourplot " << endl;
     return;
  }

  if(varindex > 4)
    varindex -= 2;


  //vector<double>  outval(nGP);
  double  outval[100];

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
          plot.triangleContourPlot(x1, x2, x3, u1, u1, u1, umin, umax, nCol);

          // contour plot for 2nd triangle
          plot.triangleContourPlot(x1, x3, x4, u1, u1, u1, umin, umax, nCol);

          count++;
      }
  }

  return;
}








void NurbsElem2DStructMixed2field::projectToKnots(bool extrapolateFlag, int vartype, int varindex, int index)
{
/*
   vals2project[0] = intVar2[indx];
   vals2project[1] = intVar2[(nGP1-1)*nivGP+indx];
   vals2project[2] = intVar2[nGP1*(nGP2-1)*nivGP+indx];
   vals2project[3] = intVar2[(nGP1*nGP2-1)*nivGP+indx];
*/

   //double outval[nGP];
   double outval[100];

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




void NurbsElem2DStructMixed2field::projectStress(int varindex, double* outval)
{

   int nivEL = nGP * nivGP;
   for(int ii=0;ii<nivEL;ii++)
     intVar2[ii] = intVar1[ii];

   int  err, isw, count, count1, ll, gp1, gp2, index;

   double  F[4], detF, F33, fact, dt, Jac, pres, cc[4][4], stre[4];

   vector<double>   N(nlbf), dN_dx(nlbf), dN_dy(nlbf);

   count  = 1;   ll = 0;  err = 0;   isw  = 3;
   dt = mpapTime.dt;

   count1 = 0;
   for(gp2=0;gp2<nGP2;gp2++)
   {
      for(gp1=0;gp1<nGP1;gp1++)
      {
          index = count1*2;

          //surf0->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, Jac);

          //surf1->deformationGradient(startindex[0], startindex[1], 1, dN_dx, dN_dy, F, detF);

          F33 = 1.0;

          matlib2d_(matDat, F, &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);

          //pres = surf2->computeValue(1, knotsAtGPs[index], knotsAtGPs[index+1]);
          
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
             cout << '\t' << "    NurbsElem2DStructMixed2field::projectStress .... : Error in 'varindex' " << endl;
             return;
          }

          count++;
          count1++;
          ll += nivGP;
     }
   }

  return;
}




void NurbsElem2DStructMixed2field::projectStrain(int vartype, int varindex, double* outval)
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






void NurbsElem2DStructMixed2field::projectIntVar(int index, double* outval)
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
void NurbsElem2DStructMixed2field::AssembleElementMatrix2(int index, MatrixXd& stiff, MatrixXd& CC)
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

void NurbsElem2DStructMixed2field::AssembleElementMatrix2(int index, MatrixXd& mat, MatrixXd& CC)
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


void NurbsElem2DStructMixed2field::AssembleElementVector2(bool firstIter, int flag, VectorXd& rhs, double* reac, VectorXd& uu)
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
      cerr << " ERROR in flag in 'NurbsElem2DStructMixed2field::AssembleElementVector2' " << endl;
      return;
   }


//  cout << " resi " << resi << endl;

  return;
}






void NurbsElem2DStructMixed2field::AssembleElementMatrix(int index, MatrixSparseArray<double>& mtx)
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





void NurbsElem2DStructMixed2field::AssembleElementVector(bool firstIter, bool flag, double* rhs, double* reac, int start1, int start2)
{
   // flag == true  ---> just external force vector
   // flag == false ---> internal load vector + contributions from nodes with specified displacement BCs

   //cout << " primvar " << endl;

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
      int aa, bb, size2, ind, *tt2;

      tt2   = &(surf2->LM[elenum2][0]);
      size2 = surf2->nsize;

      //cout << " kkkkkkkkkkkk " << endl;

      for(aa=0;aa<nsize;aa++)
      {
         if(tt[aa] != -1)
           rhs[tt[aa]] += Flocal[aa];

         // add up reaction forces
         reac[forassy[aa]] += Flocal[aa];
      }
      //cout << " qqqqqqqqqqqq " << endl;

      for(aa=0;aa<size2;aa++)
      {
         if(tt2[aa] != -1)
           rhs[tt2[aa] + start1] += resi2[aa];
      }
      
      //cout << " kkkkkkkkkkkk " << endl;

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

         for(aa=0;aa<size2;aa++)
         {
             ind = tt2[aa] + start1;

             for(bb=0;bb<nsize;bb++)
             {
                 if(tt[bb] == -1)
                 {
                   //cout << Kup(bb,aa) << '\t' << primvar[bb] << endl;
                   //cout << Kup(bb,aa) * mpapTime.dt * primvar[bb] << endl;
                   rhs[ind] -= Kup(bb,aa) * mpapTime.dt * primvar[bb];
                 }
             }
         }
      }
      tt2 = NULL;
   }

//  cout << " resi " << resi << endl;

  tt= NULL;

  return;
}


/*
void  NurbsElem2DStructMixed2field::AssembleElementMatrix(int index, Mat mtx, int start1, int start2)
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


void  NurbsElem2DStructMixed2field::AssembleElementMatrix(int index, SparseMatrixXd& mtx, int start1, int start2)
{
    int  ii, jj, nn=0, aa, bb, size2, ind, *tt1, *tt2;

    tt1 = &(surf0->LM[elenum][0]);
    tt2 = &(surf2->LM[elenum2][0]);
    size2 = surf2->nsize;
    
    /*
    cout << elenum << '\t' << elenum2 << endl;
    cout << nsize << '\t' << size2 << endl;
    cout << surf0->LM[elenum] << endl;
    cout << surf2->LM[elenum2] << endl;

   printMatrix(Klocal);
   printf("\n\n\n");
   printMatrix(Kup);
   printf("\n\n\n");
   printMatrix(Kpp);
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

        for(jj=0;jj<size2;jj++)
        {
          bb = tt2[jj];
          if(bb != -1)
          {
            bb += start1;
            mtx.coeffRef(aa, bb) += Kup(ii, jj);
            mtx.coeffRef(bb, aa) += Kpu(jj, ii);
            //mtx.coeffRef(bb, aa) += Kup(ii, jj);
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
            mtx.coeffRef(aa, bb) += Kpp(ii, jj);
          }
        }
      }
    }

  return;
}



void NurbsElem2DStructMixed2field::toPostprocess(int vartype, int varindex, int type, SparseMatrixXd&  coeffMat, VectorXd& rhsVec)
{

   double F[4], detF=0.0, F33, Jac, dt, stre[4], cc[4][4];

   int   err,  isw,  count,  count1, index, ll = 0, ii, jj, gp1, gp2, row, col;

   MatrixXd  Nlocal(nlbf,nlbf);
   VectorXd  NN(nlbf), rhslocal(nlbf), dN_dx(nlbf), dN_dy(nlbf);

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

        surf0->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), &NN(0), &dN_dx(0), &dN_dy(0), Jac);

        surf1->deformationGradient(startindex[0], startindex[1], 1, &dN_dx(0), &dN_dy(0), F, detF);

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


  return;
}







int NurbsElem2DStructMixed2field::calcError(int ind)
{
   int ii, jj, gp1, gp2, count1,  TI, index;
   
   double rad, theta, Jac, fact, dvol0;
   
   VectorXd  N(nlbf), dN_dx(nlbf), dN_dy(nlbf);
   double  v1, v2, disp1[2], disp2[2];
   
   PlateWithHole  analy(sss, matDat[1], matDat[2]);
   //ThickCylinder  analy(sss, matDat[1], matDat[2]);
   //ElasticityEx2  analy(matDat[0], matDat[1]);

   EPOINT  EP, EP2;

   double *gaussweights1 = &(surf0->gaussweights1[0]);
   double *gaussweights2 = &(surf0->gaussweights2[0]);

   double  *values1 = &(surf1->Values[0][0]);
   double  *values2 = &(surf1->Values[1][0]);
   
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

          surf0->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), &N[0], &dN_dx[0], &dN_dy[0], Jac);
      
          fact = gaussweights2[gp2] * gaussweights1[gp1] * thick * JacMultFact;

          dvol0 = Jac * fact;
          count1++;

          EP = surf0->SurfacePoint(knotsAtGPs[index], knotsAtGPs[index+1]).CalcEuclid();
          
          EP2 = surf1->SurfacePoint(knotsAtGPs[index], knotsAtGPs[index+1]).CalcEuclid();

          rad = sqrt(EP.x*EP.x + EP.y*EP.y) ;

          theta = atan2(EP.y, EP.x);
          
          //fact = analy.dispR(rad);

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
      
      VectorXd  Nbar(sizep);

      count = 1;   ll = 0;   err = 0;   isw = 3;
      dt = mpapTime.dt;

      count1 = 0;
      for(gp2=0;gp2<nGP2;gp2++)
      {
      for(gp1=0;gp1<nGP1;gp1++)
      {
          index = count1*2;

          surf0->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), &N(0), &dN_dx(0), &dN_dy(0), Jac);

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

          //stre[0] = analy.stressXX(EP.x, EP.y);
          //stre[1] = analy.stressYY(EP.x, EP.y);
          //stre[2] = 0.0;
          //stre[3] = analy.stressXY(EP.x, EP.y);

          pres = surf2->computeValueAndShanpeFns(1, knotsAtGPs[index], knotsAtGPs[index+1], &Nbar[0]);

          //
          F[0] = F[1] = F[2] = F[3] = 0.0;
          for(ii=0;ii<nlbf;ii++)
          {
            jj = tt[ii];

            bb1 = values1[jj];
            bb2 = values2[jj];

            F[0] += bb1*dN_dx[ii];
            F[2] += bb1*dN_dy[ii];
            F[1] += bb2*dN_dx[ii];
            F[3] += bb2*dN_dy[ii];
          }
          F[0] += 1.0;
          F[3] += 1.0;
          //
          
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
      double  pres, dummy, volstr;
      
      VectorXd  Nbar(sizep);

      count = 1;   ll = 0;   err = 0;   isw = 3;
      dt = mpapTime.dt;

      count1 = 0;
      for(gp2=0;gp2<nGP2;gp2++)
      {
      for(gp1=0;gp1<nGP1;gp1++)
      {
          index = count1*2;

          surf0->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), &N(0), &dN_dx(0), &dN_dy(0), Jac);
      
          fact = gaussweights2[gp2] * gaussweights1[gp1] * thick * JacMultFact;

          dvol0 = Jac * fact;

          EP = surf0->SurfacePoint(knotsAtGPs[index], knotsAtGPs[index+1]).CalcEuclid();

          rad = sqrt(EP.x*EP.x + EP.y*EP.y) ;

          theta = atan2(EP.y, EP.x);

          F1[0] = analy.strainXX(rad, theta) ;
          F1[1] = analy.strainXY(rad, theta) ;
          F1[2] = F1[1];
          F1[3] = analy.strainYY(rad, theta) ;

          //surf1->deformationGradient(startindex[0], startindex[1], 1, dN_dx, dN_dy, F, detF);
          
          //
          F2[0] = F2[1] = F2[2] = F2[3] = 0.0;
          for(ii=0;ii<nlbf;ii++)
          {
            jj = tt[ii];

            bb1 = values1[jj];
            bb2 = values2[jj];

            F2[0] += bb1*dN_dx[ii];
            F2[2] += bb1*dN_dy[ii];
            F2[1] += bb2*dN_dx[ii];
            F2[3] += bb2*dN_dy[ii];
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

        pres = surf2->computeValueAndShanpeFns(1, knotsAtGPs[index], knotsAtGPs[index+1], &Nbar(0));

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

  //if( (int) matDat[4] == 1)
    //printf(" \t element = %5d ... \t ... elemError  =   %12.6E \n " , elenum, elemError);

  return 0;
}




