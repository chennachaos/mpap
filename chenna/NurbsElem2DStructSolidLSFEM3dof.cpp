
#include <math.h>
#include "Debug.h"
#include "MpapTime.h"
#include "Plot.h"
#include "NurbsElem2DStructSolidLSFEM3dof.h"
#include "NurbsShapeFunctions.h"
#include <assert.h>
#include "ComputerTime.h"
#include "TimeFunction.h"


using namespace std;

extern ComputerTime       computerTime;
extern MpapTime mpapTime;
extern Plot plot;
extern List<TimeFunction> timeFunction;


NurbsElem2DStructSolidLSFEM3dof::NurbsElem2DStructSolidLSFEM3dof()
{
  if (debug) cout << " constructor NurbsElem2DStructSolidLSFEM3dof\n\n";
}



NurbsElem2DStructSolidLSFEM3dof::~NurbsElem2DStructSolidLSFEM3dof()
{
  if (debug) cout << " destructor NurbsElem2DStructSolidLSFEM3dof\n\n";
}



void NurbsElem2DStructSolidLSFEM3dof::createTractionDataVariable()
{
    assert(tracflag == true);
    
    if( !(tracdata.n > 0) )
    {
       tracdata.setDim(4);
       for(int ii=0;ii<4;ii++)
       {
         tracdata[ii].setDim(2);
         tracdata[ii][0] = tracdata[ii][1] = 7777.0;
       }
    }

   return;
}



int NurbsElem2DStructSolidLSFEM3dof::calcStiffnessAndResidual()
{
  if(finite)
    calcStiffnessAndResidual2();
  else
    calcStiffnessAndResidual1();

   return 0;
}


int NurbsElem2DStructSolidLSFEM3dof::calcLoadVector()
{
  if(finite)
    calcLoadVector2();
  else
    calcLoadVector1();

   return 0;
}

int NurbsElem2DStructSolidLSFEM3dof::applyDirichletBCs()
{
/*
   int  ii, jj;
   bool dispflag = false;

   for(ii=0;ii<nsize;ii++)
     stiffness_local[ii].zero();

   resi.zero();

   ListArray<VectorArray<double> >  dispdata;

   dispdata.setDim(4);
   for(ii=0;ii<4;ii++)
   {
     dispdata[ii].setDim(2);

     for(jj=0;jj<2;jj++)
       dispdata[ii][jj] = surf0->dispdata[ii][jj]; 

     printf(" \t %12.6f \t %12.6f \n", dispdata[ii][0], dispdata[ii][1]);
   }

   for(ii=0;ii<4;ii++)
   {
     dispdata[ii].setDim(2);

     for(jj=0;jj<2;jj++)
     {
       if(!CompareDoubles(dispdata[ii][jj], 7777))
       {
         dispflag = true;
         break;
       }
     }
   }

   cout << " dispflag ... = " << dispflag << endl;

   if(dispflag)
   {
      //cout << " elenum ... = " << elenum << "\n" << endl;

      double  *gausspoints1 = &(surf0->gausspoints1[0]);
      double  *gausspoints2 = &(surf0->gausspoints2[0]);
      double  *gaussweights1 = &(surf0->gaussweights1[0]);
      double  *gaussweights2 = &(surf0->gaussweights2[0]);

      double  *values1 = &(surf1->Values[0][0]);
      double  *values2 = &(surf1->Values[1][0]);
   
      int *tt = &(surf0->IEN[elenum][0]);

      double  J, Jmod, dircos[2], b1, b2, params[2], val1, val2, res;

      VectorXd f(2);
        
      int p = surf0->p, q = surf0->q, ii, jj, gp, TI, TIp1, index, nlbf2, TJ, TJp1;
      double   NN[nlbf], dN_dx[nlbf], dN_dy[nlbf], Jac, fact, fact1, fact2, ALPHA;

      ALPHA = matDat[2]; 

        // side #1
        if(!CompareDoubles(dispdata[0][0],7777))
        {
           nlbf2 = p+1;
           double  N[nlbf2];

           val1 = 0.5*uvalues[2];
           val2 = 0.5*(uvalues[0]+uvalues[1]);

           params[1] = 0.0;
           for(gp=0;gp<nGP1;gp++)   // loop over Gauss points
           {
              params[0] = val1*gausspoints1[gp] + val2;
           
              NurbsShapeFunctions2DAlg2(surf0, startindex[0], startindex[1], gausspoints1[gp], -1.0, N, J, dircos);

              surf0->ShapeFunctions(params[0], params[1], NN);

              Jmod = J * gaussweights1[gp] ;

              //printf(" res ... %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n", res, J, Jmod, Jac);

              f(0) = dispdata[0][0];
              f(1) = dispdata[0][1];
              for(ii=0;ii<nlbf;ii++)
              {
                 index = tt[ii];

                 b1 = values1[index];
                 b2 = values2[index];

                 f(0) -= b1*NN[ii];
                 f(1) -= b2*NN[ii];
              }

              printf(" f(0) ... %12.6f \t %12.6f \n", f(0), f(1));

              Jmod *= ALPHA;

              for(ii=0;ii<nlbf;ii++)
              {
                TI   = 2*ii;
                TIp1 = TI+1;

                fact1 = NN[ii]*Jmod;

                resi[TI]   += fact1*f(0);
                resi[TIp1] += fact1*f(1);

                for(jj=0;jj<nlbf;jj++)
                {
                  TJ = 2*jj;
                  TJp1 = TJ+1;

                  fact = fact1 * NN[jj] ;

                  stiffness_local[TI][TJ]      +=  fact;
                  stiffness_local[TIp1][TJp1]  +=  fact;
                }
              }
           }
           cout << " side1 done # " << endl;
        }
        // side #3
        if(!CompareDoubles(dispdata[2][0],7777))
        {
           nlbf2 = p+1;
           double   N[nlbf2];

           val1 = 0.5*uvalues[2];
           val2 = 0.5*(uvalues[0]+uvalues[1]);

           params[1] = 1.0;

           for(gp=0;gp<nGP1;gp++)   // loop over Gauss points
           {
              params[0] = val1*gausspoints1[gp] + val2;
           
              NurbsShapeFunctions2DAlg2(surf0, startindex[0], startindex[1], gausspoints1[gp], 1.0, N, J, dircos);

              surf0->ShapeFunctions(params[0], params[1], NN);

              Jmod = J * gaussweights1[gp] ;

              //printf(" res ... %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n", res, J, Jmod, Jac);
              //printf(" dircos ... %12.6f \t %12.6f \n", dircos[0], dircos[1]);

              f(0) = dispdata[2][0];
              f(1) = dispdata[2][1];
              for(ii=0;ii<nlbf;ii++)
              {
                 index = tt[ii];

                 b1 = values1[index];
                 b2 = values2[index];

                 f(0) -= b1*NN[ii];
                 f(1) -= b2*NN[ii];
              }

              printf(" f(0) ... %12.6f \t %12.6f \n", f(0), f(1));

              Jmod *= ALPHA;

              for(ii=0;ii<nlbf;ii++)
              {
                TI   = 2*ii;
                TIp1 = TI+1;

                fact1 = NN[ii]*Jmod;

                resi[TI]   += fact1*f(0);
                resi[TIp1] += fact1*f(1);

                for(jj=0;jj<nlbf;jj++)
                {
                  TJ = 2*jj;
                  TJp1 = TJ+1;

                  fact = fact1 * NN[jj] ;

                  stiffness_local[TI][TJ]      +=  fact;
                  stiffness_local[TIp1][TJp1]  +=  fact;
                }
              }
           }
           cout << " side3 done # " << endl;
        }
        // side #2
        if(!CompareDoubles(dispdata[1][0],7777) )
        {
           nlbf2 = q+1;
           double   N[nlbf2];

           val1 = 0.5*vvalues[2];
           val2 = 0.5*(vvalues[0]+vvalues[1]);

           params[0] = 1.0;

           for(gp=0;gp<nGP2;gp++)   // loop over Gauss points
           {
              params[1] = val1*gausspoints2[gp] + val2;
           
              NurbsShapeFunctions2DAlg2(surf0, startindex[0], startindex[1], 1.0, gausspoints2[gp], N, J, dircos);

              surf0->ShapeFunctions(params[0], params[1], NN);

              Jmod = J * gaussweights2[gp] ;

              //printf(" res ... %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n", res, J, Jmod, Jac);
              //printf(" dircos ... %12.6f \t %12.6f \n", dircos[0], dircos[1]);

              f(0) = dispdata[1][0];
              f(1) = dispdata[1][1];
              for(ii=0;ii<nlbf;ii++)
              {
                 index = tt[ii];

                 b1 = values1[index];
                 b2 = values2[index];

                 f(0) -= b1*NN[ii];
                 f(1) -= b2*NN[ii];
              }

              //printf(" res ... %12.6f \n", res);

              Jmod *= ALPHA;

              for(ii=0;ii<nlbf;ii++)
              {
                TI   = 2*ii;
                TIp1 = TI+1;

                fact1 = NN[ii]*Jmod;

                resi[TI]   += fact1*f(0);
                resi[TIp1] += fact1*f(1);

                for(jj=0;jj<nlbf;jj++)
                {
                  TJ = 2*jj;
                  TJp1 = TJ+1;

                  fact = fact1 * NN[jj] ;

                  stiffness_local[TI][TJ]      +=  fact;
                  stiffness_local[TIp1][TJp1]  +=  fact;
                }
              }
           }
           cout << " side2 done # " << endl;
        }
        // side #4
        if(!CompareDoubles(dispdata[3][0],7777) )
        {
           nlbf2 = q+1;
           double   N[nlbf2];

           val1 = 0.5*vvalues[2];
           val2 = 0.5*(vvalues[0]+vvalues[1]);

           params[0] = 0.0;

           for(gp=0;gp<nGP2;gp++)   // loop over Gauss points
           {
              params[1] = val1*gausspoints2[gp] + val2;
           
              NurbsShapeFunctions2DAlg2(surf0, startindex[0], startindex[1], -1.0, gausspoints2[gp], N, J, dircos);

              surf0->ShapeFunctions(params[0], params[1], NN);

              Jmod = J * gaussweights2[gp] ;

              //printf(" res ... %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n", res, J, Jmod, Jac);
              //printf(" dircos ... %12.6f \t %12.6f \n", dircos[0], dircos[1]);

              f(0) = dispdata[3][0];
              f(1) = dispdata[3][1];
              for(ii=0;ii<nlbf;ii++)
              {
                 index = tt[ii];

                 b1 = values1[index];
                 b2 = values2[index];

                 f(0) -= b1*NN[ii];
                 f(1) -= b2*NN[ii];
              }

              //printf(" res ... %12.6f \n", res);

              Jmod *= ALPHA;

              for(ii=0;ii<nlbf;ii++)
              {
                TI   = 2*ii;
                TIp1 = TI+1;

                fact1 = NN[ii]*Jmod;

                resi[TI]   += fact1*f(0);
                resi[TIp1] += fact1*f(1);

                for(jj=0;jj<nlbf;jj++)
                {
                  TJ = 2*jj;
                  TJp1 = TJ+1;

                  fact = fact1 * NN[jj] ;

                  stiffness_local[TI][TJ]      +=  fact;
                  stiffness_local[TIp1][TJp1]  +=  fact;
                }
              }
           }
           cout << " side4 done # " << endl;
        }
//printForceVector();
//printf("\n\n");
//printStiffnessMatrix();
    }
*/
  return 0;
}




/*
int NurbsElem2DStructSolidLSFEM3dof::calcStiffnessAndResidual1()
{
   int ii, jj, gp1, gp2, TI, TIp1, TIp2, count1, index;

   double  F[4], detF=0.0, F33, fact, dvol, dvol0, Jac, dt, ff, lambda, BULK, mu, pres, nu;
   
   double  b1, b2, b3, trgradu, ci;

   double  N[nlbf], dN_dx[nlbf], d2N_dx2[nlbf], dN_dy[nlbf], d2N_dy2[nlbf];
   VectorXd  f(3), dp(3), Du(3);
   MatrixXd  D(nsize,3);
   
   BULK = matDat[0];
   mu = matDat[1];
   lambda = BULK - 2.0*mu/3.0;

   BULK = lambda+mu;

   ff = 1.0;

   for(ii=0;ii<nsize;ii++)
     stiffness_local[ii].zero();
   
   resi.zero();

   double *gaussweights = &(surf0->gaussweights[0]);

   double  *values1 = &(surf1->Values[0][0]);
   double  *values2 = &(surf1->Values[1][0]);
   double  *values3 = &(surf1->Values[2][0]);
   
   int *tt = &(surf0->IEN[elenum][0]);

   count1 = 0;
   for(gp2=0;gp2<nGP2;gp2++)
   {
   for(gp1=0;gp1<nGP1;gp1++)
   {
        index = count1*2;

        surf0->ShapeFunDerivatives2(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, d2N_dx2, d2N_dy2, Jac);

        dvol0 = Jac * gaussweights[count1] * JacMultFact;
        count1++;

        //surf1->deformationGradient(startindex[0], startindex[1], 1, &dN_dx(0), &dN_dy(0), F, detF);

        //for(ii=0;ii<nlbf;ii++)
          //printf(" \t %14.8f\t %14.8f\t%14.8f\t %14.8f\t%14.8f\n", N[ii], dN_dx[ii], dN_dy[ii], d2N_dx2[ii], d2N_dy2[ii]);
        //printf("\n\n");

        f.setZero();
        dp.setZero();
        Du.setZero();
        trgradu = pres = 0.0;

        f(0) = -bforce[0] * rho0 ;
        f(1) = -bforce[1] * rho0 ;

        for(ii=0;ii<nlbf;ii++)
        {
           index = tt[ii];
             
           b1 = values1[index];
           b2 = values2[index];
           b3 = values3[index];
           
           fact = d2N_dx2[ii]+d2N_dy2[ii];

           Du(0)  += (b1*fact);
           Du(1)  += (b2*fact);
           
           trgradu  +=  (b1*dN_dx[ii]+b2*dN_dy[ii]);
           
           pres += (b3 * N[ii]);
           
           dp(0)  += (b3 * dN_dx[ii]);
           dp(1)  += (b3 * dN_dy[ii]);
           
           /////////////////////////////////////

           TI   =  3*ii;
           TIp1 =  TI+1;
           TIp2 =  TI+2;

           b1 = dN_dx[ii];
           b2 = dN_dy[ii];

           ci = mu * fact;

           D(TI,0)   = ci;
           D(TIp1,0) = 0.0;
           D(TIp2,0) = ff*b1;

           D(TI,1)   = 0.0;
           D(TIp1,1) = ci;
           D(TIp2,1) = ff*b2;

           D(TI,2)   =  b1;
           D(TIp1,2) =  b2;
           D(TIp2,2) = -N[ii]/BULK;
        }

        //printf(" %14.8f \t%14.8f \t%14.8f \t%14.8f \t%14.8f \t%14.8f \t%14.8f \t%14.8f\n\n", ff, Du(0), Du(1), dp(0), dp(1), f(0), f(1), pres);

        f -= (mu*Du + ff*dp);
        f(2) = -(trgradu - pres/BULK);

        for(ii=0;ii<nsize;ii++)
        {
           resi[ii] += dvol0*(D(ii,0)*f(0)+D(ii,1)*f(1)+D(ii,2)*f(2));

           for(jj=0;jj<nsize;jj++)
             stiffness_local[ii][jj]  +=  dvol0*(D(ii,0)*D(jj,0)+D(ii,1)*D(jj,1) + D(ii,2)*D(jj,2)) ;
        }
    }//gp1
    }//gp2

//printStiffnessMatrix();
//printf("\n\n");
//printForceVector();

   return 0;
}
*/


int NurbsElem2DStructSolidLSFEM3dof::toComputeInfSupCondition()
{
/*
   // to compute the inf-sup constant
   
   int  err, isw, count, count1, ll, ii, jj, twoI, twoJ, twoIp1, twoJp1, index, gp1, gp2, mm;
   int  sizep = surf2->nlbf;

   double  F[4], detF, F33, fact, dvol0, dt, Jac, dummy, pres, b1, b2, b3, b4, b5, volstr;
   double  cc[4][4], stre[4], bc[2][3], Idev[4][4], cctmp[4][4];
   double  N[nlbf], dN_dx[nlbf], dN_dy[nlbf], Nbar[sizep];
   double  BULK = matDat[0];
   double  eta = elmDat[7];

   Idev2D(Idev);

   double *gaussweights = &(surf0->gaussweights[0]);

   for(ii=0;ii<nsize;ii++)
     stiffness_local[ii].zero();

   resi.zero();
   resi2.zero();

   Kup.resize(nsize, nlbf);
   Kup.setZero();

   Kpp.resize(nlbf, nlbf);
   Kpp.setZero();

   count = 1;   ll = 0;   err = 0;   isw = 3;
   dt = mpapTime.dt;

   count1 = 0;
   for(gp2=0;gp2<nGP2;gp2++)
   {
      for(gp1=0;gp1<nGP1;gp1++)
      {
          index = count1*2;

          surf0->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, Jac);
          
          surf0->ShapeFunctions(knotsAtGPs[index], knotsAtGPs[index+1], N);

          dvol0 = Jac * gaussweights[count1] * JacMultFact;

          //==============================================
          // CALCULATE TANGENT STIFFNESS and RESIDUAL
          //==============================================

          for(ii=0;ii<nlbf;ii++)
          {
             twoI   = 2*ii;
             twoIp1 = twoI+1;

             b1 = dN_dx[ii]*dvol0;
             b2 = dN_dy[ii]*dvol0;
             b3 = N[ii]*dvol0;

             for(jj=0;jj<nlbf;jj++)
             {
                twoJ   = 2*jj;
                twoJp1 = twoJ+1;

                fact = b1*dN_dx[jj] + b2*dN_dy[jj] + eta*b3*N[jj];

                stiffness_local[twoI][twoJ]     += fact;
                stiffness_local[twoI][twoJp1]   += 0.0;
                stiffness_local[twoIp1][twoJ]   += 0.0;
                stiffness_local[twoIp1][twoJp1] += fact;
             }
             // compute Kup and Kpp matrices

             for(jj=0;jj<nlbf;jj++)
             {
                Kup(twoI, jj)   += ( b1 * N[jj] );
                Kup(twoIp1, jj) += ( b2 * N[jj] );

                Kpp(ii,jj) += ( b3 * N[jj] );
             }
          }

          count1++;
     }
   }

//  printStiffnessMatrix();
//  printf("\n\n");
//  printForceVector();
*/
  return 0;
}


//
int NurbsElem2DStructSolidLSFEM3dof::calcStiffnessAndResidual1()
{
/*
   int ii, jj, kk, mm, ll, gp1, gp2, TI, TIp1, TIp2, TJ, TJp1, TJp2, count1, index, isw, count, err;

   double  F[4], detF=0.0, F33, fact, dvol, dvol0, Jac, dt, ff, lambda, BULK, mu, pres, nu;
   
   double  b1, b2, b3, c1, c2, c3, trgradu, ci, volstr, dummy;

   double  cc[4][4], strain[4], stre[4], bc[2][3], Idev[4][4], cctmp[4][4];
   double  N[nlbf], dN_dx[nlbf], d2N_dx2[nlbf], dN_dy[nlbf], d2N_dy2[nlbf], d2N_dxy[nlbf], d2N_dyx[nlbf];
   VectorXd  f(3), dp(3), Du(3);
   MatrixXd  D(nsize,3);
   
   Idev2D(Idev);
   
   b1 =  4.0/3.0;
   b2 = -2.0/3.0;
   
   cctmp[0][0] = b1;  cctmp[0][1] = b2;  cctmp[0][2] = b2;  cctmp[0][3] = 0.0;
   cctmp[1][0] = b2;  cctmp[1][1] = b1;  cctmp[1][2] = b2;  cctmp[1][3] = 0.0;
   cctmp[2][0] = b2;  cctmp[2][1] = b2;  cctmp[2][2] = b1;  cctmp[2][3] = 0.0;
   cctmp[3][0] = 0.0; cctmp[3][1] = 0.0; cctmp[3][2] = 0.0; cctmp[3][3] = 1.0;


   BULK = matDat[0];
   mu = matDat[1];
   lambda = BULK - 2.0*mu/3.0;

   for(ii=0;ii<nsize;ii++)
     stiffness_local[ii].zero();
   
   resi.zero();

   double *gaussweights = &(surf0->gaussweights[0]);

   double  *values1 = &(surf1->Values[0][0]);
   double  *values2 = &(surf1->Values[1][0]);
   double  *values3 = &(surf1->Values[2][0]);
   
   int *tt = &(surf0->IEN[elenum][0]);

   count = 1;   ll = 0;   err = 0;   isw = 3;
   dt = mpapTime.dt;

   bforce[0] = elmDat[5] ;
   bforce[1] = elmDat[6] ;
   rho0      = elmDat[7] ;


   count1 = 0;
   for(gp2=0;gp2<nGP2;gp2++)
   {
   for(gp1=0;gp1<nGP1;gp1++)
   {
        index = count1*2;

        surf0->ShapeFunDerivatives2(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, d2N_dx2, d2N_dy2, d2N_dxy, d2N_dyx, Jac);

        dvol0 = Jac * gaussweights[count1] * JacMultFact;
        count1++;

        surf1->deformationGradient(startindex[0], startindex[1], 1, dN_dx, dN_dy, F, detF);

        //for(ii=0;ii<nlbf;ii++)
          //printf(" \t %14.8f\t %14.8f\t%14.8f\t %14.8f\t%14.8f\n", N[ii], dN_dx[ii], dN_dy[ii], d2N_dx2[ii], d2N_dy2[ii]);
        //printf("\n\n");

        trgradu = pres = 0.0;

        //F[1] = F[2] = 0.0;
        //F[0] = F[3] = 1.0;
        for(ii=0;ii<nlbf;ii++)
        {
           index = tt[ii];
             
           b1 = values1[index];
           b2 = values2[index];
           b3 = values3[index];
           
           //F[0] += dN_dx[ii]*b1;
           //F[1] += dN_dy[ii]*b1;
           //F[2] += dN_dx[ii]*b2;
           //F[3] += dN_dy[ii]*b2;

           pres += N[ii]*b3;
        }

        //trgradu  = F[0] + F[3];
        //detF = F[0]*F[3] - F[1]*F[2];

          F33 = 1.0;

          matlib2d_(matDat, F, &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);
          if(err !=0)           return 1;

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

          strain[0] = F[0] - 1.0;
          strain[1] = F[3] - 1.0;
          strain[2] = F33 - 1.0;
          strain[3] = 0.5*(F[1]+F[2]);
          
          volstr = strain[0] + strain[1] + strain[2];
          volstr = F[0] + F[3] - 2.0;
          
          //stre[0] = 2.0*mu*strain[0] + lambda*volstr;
          //stre[1] = 2.0*mu*strain[1] + lambda*volstr;
          //stre[2] = 2.0*mu*strain[2] + lambda*volstr;
          //stre[3] = 2.0*mu*strain[3];
          
          dummy = pres - (stre[0]+stre[1]+stre[2])/3.0 ;

          stre[0] += dummy;
          stre[1] += dummy;
          stre[2] += dummy;

          //for(ii=0;ii<4;ii++)
            //stre[ii] *= 2.0;
          //for(ii=0;ii<4;ii++)
          //{
            //for(jj=0;jj<4;jj++)
              //cc[ii][jj] = mu*cctmp[ii][jj];
          //}

          //printf(" stresses ");        printf("\t%12.8f\t%12.8f\t%12.8f\t%12.8f\n\n", stre[0], stre[1], stre[2], pres);

          //==============================================
          // CALCULATE TANGENT STIFFNESS and RESIDUAL
          //==============================================
          
          
          f(0) =  bforce[0]*rho0;
          f(1) =  bforce[1]*rho0;
          f(2) =  0.0;

          for(ii=0;ii<nlbf;ii++)
          {
             b1 = dN_dx[ii]*dvol0;
             b2 = dN_dy[ii]*dvol0;
             b3 = N[ii] * dvol0;

             bc[0][0] = (b1 * cc[0][0] + b2 * cc[3][0]);
             bc[0][1] = (b1 * cc[0][1] + b2 * cc[3][1]);
             bc[0][2] = (b1 * cc[0][3] + b2 * cc[3][3]);

             bc[1][0] = (b2 * cc[1][0] + b1 * cc[3][0]);
             bc[1][1] = (b2 * cc[1][1] + b1 * cc[3][1]);
             bc[1][2] = (b2 * cc[1][3] + b1 * cc[3][3]);

             TI   = 3*ii;
             TIp1 = TI+1;
             TIp2 = TI+2;

             resi[TI]   += (b3*f(0) - b1*stre[0] - b2*stre[3]) ;
             resi[TIp1] += (b3*f(1) - b1*stre[3] - b2*stre[1]) ;
             resi[TIp2] -= (b3*(volstr - pres/BULK));

             for(jj=0;jj<nlbf;jj++)
             {
                TJ   = 3*jj;
                TJp1 = TJ+1;
                TJp2 = TJ+2;

                c1 = dN_dx[jj];
                c2 = dN_dy[jj];
                c3 = N[jj];

                stiffness_local[TI][TJ]     +=  (bc[0][0] * c1 + bc[0][2] * c2) ;
                stiffness_local[TI][TJp1]   +=  (bc[0][1] * c2 + bc[0][2] * c1) ;
                stiffness_local[TIp1][TJ]   +=  (bc[1][0] * c1 + bc[1][2] * c2) ;
                stiffness_local[TIp1][TJp1] +=  (bc[1][1] * c2 + bc[1][2] * c1) ;

                // pressure term
                stiffness_local[TI][TJp2]   += (b1*c3);
                stiffness_local[TIp1][TJp2] += (b2*c3);

                // constraint equation
                stiffness_local[TIp2][TJ]   += (b3*c1);
                stiffness_local[TIp2][TJp1] += (b3*c2);
                //stiffness_local[TIp2][TJp2] += 0.0;
                stiffness_local[TIp2][TJp2] -= b3*c3/BULK;
             }
          }

          //count++;
          //count1++;
          //ll += nivGP;

    }//gp1
    }//gp2

  //printStiffnessMatrix();
  //printf("\n\n");
  //printForceVector();
*/
  return 0;
}
//


/*
int NurbsElem2DStructSolidLSFEM3dof::calcLoadVector1()
{
   if(tracflag)
   {
      for(int ii=0;ii<nsize;ii++)
        stiffness_local[ii].zero();

      resi.zero();

      double *gausspoints1 = &(surf0->gausspoints1[0]);
      double *gausspoints2 = &(surf0->gausspoints2[0]);
      double *gaussweights1 = &(surf0->gaussweights1[0]);
      double *gaussweights2 = &(surf0->gaussweights2[0]);
    
      VectorXd  trac(2), res(2);
      MatrixXd  D(nsize, 2), F(2,2);
    
      double  *values1 = &(surf1->Values[0][0]);
      double  *values2 = &(surf1->Values[1][0]);
      double  *values3 = &(surf1->Values[2][0]);
   
      int *tt = &(surf0->IEN[elenum][0]);

      double  J, J1, Jmod, dircos[2], BULK, mu, a1, a2, b1, b2, b3, params[2], val1, val2, d, pres, lambda;
        
      BULK = matDat[0];
      mu = matDat[1];

      lambda = BULK - 2.0*mu/3.0;

      BULK = lambda+mu;

      d = 2.0;

      a1 = 1.0/d;
      a2 = 1.0- a1;
        
      int p = surf0->p, q = surf0->q, ii, jj, gp, TI, TIp1, TIp2, index, nlbf2, TJ, TJp1, TJp2;
      double   NN[nlbf], dN_dx[nlbf], dN_dy[nlbf], Jac, fact, fact1, fact2, ALPHA, ALPHA1;

      ALPHA = matDat[2];
      ALPHA1 = matDat[3];

      ALPHA1 = ALPHA = 1.0;

      bool pout = false;

      pout = (bool) matDat[3];

      // side #1
      if(!CompareDoubles(tracdata[0][0],7777) || !CompareDoubles(tracdata[0][1],7777))
      {
          nlbf2 = p+1;
          double  N[nlbf2];

          val1 = 0.5*uvalues[2];
          val2 = 0.5*(uvalues[0]+uvalues[1]);

          params[1] = 0.0;
          for(gp=0;gp<nGP1;gp++)   // loop over Gauss points
          {
              params[0] = val1*gausspoints1[gp] + val2;
           
              NurbsShapeFunctions2DAlg2(surf0, startindex[0], startindex[1], gausspoints1[gp], -1.0, N, J, dircos);

              surf0->ShapeFunDerivatives3(&(startindex[0]), params, NN, dN_dx, dN_dy, Jac);

              Jmod = J * gaussweights1[gp];

              //dircos[0] *= -1.0;
              //dircos[1] *= -1.0;

              res(0) = tracdata[0][0] * (-dircos[0]) + tracdata[0][1] * (-dircos[1]);
              res(1) = tracdata[0][0] * (-dircos[1]) + tracdata[0][1] * (dircos[0]);

              if(pout)
              {
                printf(" tracX and tracY ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n", res(0), res(1), J, Jmod, Jac);
                printf(" dircos ... %12.6f \t %12.6f \n", dircos[0], dircos[1]);
              }
              
              F.setZero();
              pres = 0.0;
              for(ii=0;ii<nlbf;ii++)
              {
                 index = tt[ii];

                 b1 = values1[index];
                 b2 = values2[index];
                 b3 = values3[index];

                 F(0,0) += b1 * dN_dx[ii];
                 F(0,1) += b1 * dN_dy[ii];
                 F(1,0) += b2 * dN_dx[ii];
                 F(1,1) += b2 * dN_dy[ii];
                 
                 pres += b3 * NN[ii];

                 TI = 3 * ii;
                 TIp1 = TI+1;
                 TIp2 = TI+2;

                 D(TI,0)   = mu*( 2.0*a2*dN_dx[ii]*dircos[0] + dN_dy[ii]*dircos[1]);
                 D(TIp1,0) = mu*(-2.0*a1*dN_dy[ii]*dircos[0] + dN_dx[ii]*dircos[1]);
                 D(TIp2,0) = NN[ii]*dircos[0];

                 D(TI,1)   = mu*(dN_dy[ii]*dircos[0] - 2.0*a1*dN_dx[ii]*dircos[1]);
                 D(TIp1,1) = mu*(dN_dx[ii]*dircos[0] + 2.0*a2*dN_dy[ii]*dircos[1]);
                 D(TIp2,1) = NN[ii]*dircos[1];
              }
              //printf("\n\n");
              
              res(0) -= ((2.0*mu*(a2*F(0,0) - a1*F(1,1)) + pres)*dircos[0] + mu*(F(0,1)+F(1,0))*dircos[1]);
              res(1) -= (mu*(F(0,1)+F(1,0))*dircos[0] + (2.0*mu*(a2*F(1,1) - a1*F(0,0))+pres)*dircos[1]);

              Jmod *= ALPHA1;
              for(ii=0;ii<nsize;ii++)
              {
                resi[ii] += Jmod*(D(ii,0)*res(0)+D(ii,1)*res(1));

                for(jj=0;jj<nsize;jj++)
                  stiffness_local[ii][jj]  +=  Jmod*(D(ii,0)*D(jj,0)+D(ii,1)*D(jj,1)) ;
              }
          }
          if(pout) cout << " side1 done " << endl;
        }
        // side #3
        if(!CompareDoubles(tracdata[2][0],7777) || !CompareDoubles(tracdata[2][1],7777))
        {
           nlbf2 = p+1;
           double   N[nlbf2];

           val1 = 0.5*uvalues[2];
           val2 = 0.5*(uvalues[0]+uvalues[1]);

           params[1] = 1.0;

           for(gp=0;gp<nGP1;gp++)   // loop over Gauss points
           {
              params[0] = val1*gausspoints1[gp] + val2;
           
              NurbsShapeFunctions2DAlg2(surf0, startindex[0], startindex[1], gausspoints1[gp], 1.0, N, J, dircos);

              surf0->ShapeFunDerivatives3(&(startindex[0]), params, NN, dN_dx, dN_dy, Jac);

              Jmod = J * gaussweights1[gp] ;

              dircos[0] *= -1.0;
              dircos[1] *= -1.0;
              
              res(0) = tracdata[2][0] * (-dircos[0]) + tracdata[2][1] * (-dircos[1]);
              res(1) = tracdata[2][0] * (-dircos[1]) + tracdata[2][1] * (dircos[0]);

              if(pout)
              {
                printf(" tracX and tracY ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n", res(0), res(1), J, Jmod, Jac);
                printf(" dircos ... %12.6f \t %12.6f \n", dircos[0], dircos[1]);
              }
              
              F.setZero();
              pres = 0.0;
              for(ii=0;ii<nlbf;ii++)
              {
                 index = tt[ii];

                 b1 = values1[index];
                 b2 = values2[index];
                 b3 = values3[index];

                 F(0,0) += b1 * dN_dx[ii];
                 F(0,1) += b1 * dN_dy[ii];
                 F(1,0) += b2 * dN_dx[ii];
                 F(1,1) += b2 * dN_dy[ii];
                 
                 pres += b3 * NN[ii];

                 TI = 3 * ii;
                 TIp1 = TI+1;
                 TIp2 = TI+2;

                 D(TI,0)   = mu*( 2.0*a2*dN_dx[ii]*dircos[0] + dN_dy[ii]*dircos[1]);
                 D(TIp1,0) = mu*(-2.0*a1*dN_dy[ii]*dircos[0] + dN_dx[ii]*dircos[1]);
                 D(TIp2,0) = NN[ii]*dircos[0];

                 D(TI,1)   = mu*(dN_dy[ii]*dircos[0] - 2.0*a1*dN_dx[ii]*dircos[1]);
                 D(TIp1,1) = mu*(dN_dx[ii]*dircos[0] + 2.0*a2*dN_dy[ii]*dircos[1]);
                 D(TIp2,1) = NN[ii]*dircos[1];
              }
              
              res(0) -= ((2.0*mu*(a2*F(0,0) - a1*F(1,1)) + pres)*dircos[0] + mu*(F(0,1)+F(1,0))*dircos[1]);
              res(1) -= (mu*(F(0,1)+F(1,0))*dircos[0] + (2.0*mu*(a2*F(1,1) - a1*F(0,0))+pres)*dircos[1]);

              Jmod *= ALPHA1;
              for(ii=0;ii<nsize;ii++)
              {
                resi[ii] += Jmod*(D(ii,0)*res(0)+D(ii,1)*res(1));

                for(jj=0;jj<nsize;jj++)
                  stiffness_local[ii][jj]  +=  Jmod*(D(ii,0)*D(jj,0)+D(ii,1)*D(jj,1)) ;
              }
           }
           if(pout) cout << " side3 done " << endl;
        }
        // side #2
        if(!CompareDoubles(tracdata[1][0],7777) || !CompareDoubles(tracdata[1][1],7777))
        {
           nlbf2 = q+1;
           double   N[nlbf2];

           val1 = 0.5*vvalues[2];
           val2 = 0.5*(vvalues[0]+vvalues[1]);

           params[0] = 1.0;

           for(gp=0;gp<nGP2;gp++)   // loop over Gauss points
           {
              params[1] = val1*gausspoints2[gp] + val2;
           
              NurbsShapeFunctions2DAlg2(surf0, startindex[0], startindex[1], 1.0, gausspoints2[gp], N, J, dircos);

              surf0->ShapeFunDerivatives3(&(startindex[0]), params, NN, dN_dx, dN_dy, Jac);

              Jmod = J * gaussweights2[gp] ;
              
              res(0) = tracdata[1][0] * (-dircos[0]) + tracdata[1][1] * (-dircos[1]);
              res(1) = tracdata[1][0] * (-dircos[1]) + tracdata[1][1] * (dircos[0]);

              //res(0) = tracdata[1][0];
              //res(1) = tracdata[1][1];

              if(pout)
              {
                printf(" tracX and tracY ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n", res(0), res(1), J, Jmod, Jac);
                printf(" dircos ... %12.6f \t %12.6f \n", dircos[0], dircos[1]);
              }
              
              F.setZero();
              pres = 0.0;
              for(ii=0;ii<nlbf;ii++)
              {
                 index = tt[ii];

                 b1 = values1[index];
                 b2 = values2[index];
                 b3 = values3[index];

                 F(0,0) += b1 * dN_dx[ii];
                 F(0,1) += b1 * dN_dy[ii];
                 F(1,0) += b2 * dN_dx[ii];
                 F(1,1) += b2 * dN_dy[ii];
                 
                 pres += b3 * NN[ii];

                 TI = 3 * ii;
                 TIp1 = TI+1;
                 TIp2 = TI+2;

                 D(TI,0)   = mu*( 2.0*a2*dN_dx[ii]*dircos[0] + dN_dy[ii]*dircos[1]);
                 D(TIp1,0) = mu*(-2.0*a1*dN_dy[ii]*dircos[0] + dN_dx[ii]*dircos[1]);
                 D(TIp2,0) = NN[ii]*dircos[0];

                 D(TI,1)   = mu*(dN_dy[ii]*dircos[0] - 2.0*a1*dN_dx[ii]*dircos[1]);
                 D(TIp1,1) = mu*(dN_dx[ii]*dircos[0] + 2.0*a2*dN_dy[ii]*dircos[1]);
                 D(TIp2,1) = NN[ii]*dircos[1];
              }
              
              res(0) -= ((2.0*mu*(a2*F(0,0) - a1*F(1,1)) + pres)*dircos[0] + mu*(F(0,1)+F(1,0))*dircos[1]);
              res(1) -= (mu*(F(0,1)+F(1,0))*dircos[0] + (2.0*mu*(a2*F(1,1) - a1*F(0,0))+pres)*dircos[1]);

              Jmod *= ALPHA1;
              for(ii=0;ii<nsize;ii++)
              {
                resi[ii] += Jmod*(D(ii,0)*res(0)+D(ii,1)*res(1));

                for(jj=0;jj<nsize;jj++)
                  stiffness_local[ii][jj]  +=  Jmod*(D(ii,0)*D(jj,0)+D(ii,1)*D(jj,1)) ;
              }
           }
           if(pout) cout << " side2 done " << endl;
        }
        // side #4
        if(!CompareDoubles(tracdata[3][0],7777) || !CompareDoubles(tracdata[3][1],7777))
        {
           nlbf2 = q+1;
           double   N[nlbf2];

           val1 = 0.5*vvalues[2];
           val2 = 0.5*(vvalues[0]+vvalues[1]);

           params[0] = 0.0;

           for(gp=0;gp<nGP2;gp++)   // loop over Gauss points
           {
              params[1] = val1*gausspoints2[gp] + val2;
           
              NurbsShapeFunctions2DAlg2(surf0, startindex[0], startindex[1], -1.0, gausspoints2[gp], N, J, dircos);

              surf0->ShapeFunDerivatives3(&(startindex[0]), params, NN, dN_dx, dN_dy, Jac);

              Jmod = J * gaussweights2[gp] ;

              dircos[0] *= -1.0;
              dircos[1] *= -1.0;

              res(0) = tracdata[3][0] * (-dircos[0]) + tracdata[3][1] * (-dircos[1]);
              res(1) = tracdata[3][0] * (-dircos[1]) + tracdata[3][1] * (dircos[0]);

              //printf(" tracX and tracY ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n", res(0), res(1), J, Jmod, Jac);
              //printf(" dircos ... %12.6f \t %12.6f \n", dircos[0], dircos[1]);

              res.setZero();
              pres = 0.0;
              for(ii=0;ii<nlbf;ii++)
              {
                 index = tt[ii];

                 b1 = values1[index];
                 b2 = values2[index];
                 b3 = values3[index];

                 res(0) -= b1*NN[ii];
                 res(1) -= b2*NN[ii];
                 pres -= b3*NN[ii];
              }

              Jmod *= ALPHA;

              for(ii=0;ii<nlbf;ii++)
              {
                TI   = 3*ii;
                TIp1 = TI+1;
                TIp2 = TI+2;

                fact1 = NN[ii]*Jmod;

                resi[TI]   += fact1*res(0);
                resi[TIp1] += fact1*res(1);
                //resi[TIp2] += fact1*pres;

                for(jj=0;jj<nlbf;jj++)
                {
                  TJ   = 3*jj;
                  TJp1 = TJ+1;
                  TJp2 = TJ+2;

                  fact = fact1 * NN[jj];

                  stiffness_local[TI][TJ]      +=  fact;
                  stiffness_local[TIp1][TJp1]  +=  fact;
                  //stiffness_local[TIp2][TJp2]  +=  fact;
                }
              }           
           }
           if(pout) cout << " side4 done " << endl;
        }
if(pout) printForceVector();
    }
//printStiffnessMatrix();
//printf("\n\n");
//printForceVector();

  return 0;
}
*/



//
int NurbsElem2DStructSolidLSFEM3dof::calcLoadVector1()
{
/*
  //cout << " AAAAAAAAAA " << endl;

    // initialize resi vector
    resi.zero();

   double *gausspoints1 = &(surf0->gausspoints1[0]);
   double *gausspoints2 = &(surf0->gausspoints2[0]);
   double *gaussweights1 = &(surf0->gaussweights1[0]);
   double *gaussweights2 = &(surf0->gaussweights2[0]);

    //   1.) BODY FORCES
    //============================================

    if(!CompareDoubles(bforce[0],0.0) || !CompareDoubles(bforce[1],0.0))
    {
       double  NN[nlbf], J=0.0, dvol0, dvol, fact1, fact2, fact3;
       int  count=0, index, gp2, gp1, ii;

       // loop over Gauss points
       for(gp2=0;gp2<nGP2;gp2++)
       {
       for(gp1=0;gp1<nGP1;gp1++)
       {
              index = 2*count;
              NurbsShapeFunctions2DAlg3(surf0, startindex[0], startindex[1], knotsAtGPs[index], knotsAtGPs[index+1], NN, J);
              count++;

              dvol0 = J * gaussweights1[gp1] * gaussweights2[gp2] * thick;

              //   for axisymmetric problems compute radius
              if(axsy)
              {
                 double rad0 = 0.0;
                 int loc_num=0;
                 for(int jj=0; jj<=surf0->q; jj++)
                 {
                    for(ii=0; ii<=surf0->p; ii++)
                    {
                       rad0 += (NN[loc_num] * surf0->Pw[startindex[0]+ii][startindex[1]+jj].CalcEuclid().x);
                       loc_num++;
                    }
                 }
                 dvol0 *= (twoPI * rad0);
              }

              fact3 = rho0 * dvol0 ;
              fact1 = bforce[0] * fact3 ;
              fact2 = bforce[1] * fact3 ;

              for(ii=0;ii<nlbf;ii++)
              {
                 resi[3*ii]   +=  NN[ii] * fact1 ;
                 resi[3*ii+1] +=  NN[ii] * fact2 ;
              }
       }//gp1
       }//gp2
    }

     //   2.) FORCES DUE TO SURFACE TRACTION
     //============================================
     // normal traction into the surface is taken as positive

     
     if(tracflag)
     {
        //
        cout << "       elem... : " << elenum << endl;
        cout << tracdata[0][0] << '\t' << tracdata[0][1] << '\t' << CompareDoubles(tracdata[0][0],0.0)<< '\t' << CompareDoubles(tracdata[0][1],0.0) << endl;
        cout << tracdata[1][0] << '\t' << tracdata[1][1] << '\t' << CompareDoubles(tracdata[1][0],0.0)<< '\t' << CompareDoubles(tracdata[1][1],0.0) << endl;
        cout << tracdata[2][0] << '\t' << tracdata[2][1] << '\t' << CompareDoubles(tracdata[2][0],0.0)<< '\t' << CompareDoubles(tracdata[2][1],0.0) << endl;
        cout << tracdata[3][0] << '\t' << tracdata[3][1] << '\t' << CompareDoubles(tracdata[3][0],0.0)<< '\t' << CompareDoubles(tracdata[3][1],0.0) << endl;
        //

        double dvol0=0.0, fact1, J, Jmod, rad=0.0, dircos[2], tracX, tracY, param;
        int index=0, ngbf1, ngbf2;

        int p = surf0->p, q = surf0->q, ind1, ind2, ii, gp;

        ListArray<CPOINT>  Pw1;
        ListArray<EPOINT>  SKL;

        EPOINT  EP1, EP2, Normal;

        ngbf1 = surf0->ngbf1;
        ngbf2 = surf0->ngbf2;


        // side #1
        if(!CompareDoubles(tracdata[0][0],7777.0) || !CompareDoubles(tracdata[0][1],7777.0))
        {
            double NN[p+1];

            Pw1.setDim(ngbf1);

            for(ii=0;ii<ngbf1;ii++)
               Pw1[ii] = surf0->Pw[ii][0];

           NurbsCURVE   curve_temp(Pw1, surf0->U, p);

           for(gp=0;gp<nGP1;gp++)   // loop over Gauss points
           {
              if( finite && followerLoadFlag )
                 NurbsShapeFunctions2DAlg2(surf1, startindex[0], startindex[1], gausspoints1[gp], -1.0, NN, J, dircos);
              else
                 NurbsShapeFunctions2DAlg2(surf0, startindex[0], startindex[1], gausspoints1[gp], -1.0, NN, J, dircos);

              Jmod = J * gaussweights1[gp] * thick;

              //   for axisymmetric problems compute radius
              if(axsy)
              {
                 if(finite && followerLoadFlag )
                    rad  = surf1->SurfacePoint(knotsAtGPs[2*gp], vvalues[0]).CalcEuclid().x;
                 else
                    rad  = surf0->SurfacePoint(knotsAtGPs[2*gp], vvalues[0]).CalcEuclid().x;

                 Jmod *= (twoPI * rad);
              }

              tracX = tracdata[0][0] * (-dircos[0]) + tracdata[0][1] * (-dircos[1]);
              tracY = tracdata[0][0] * (-dircos[1]) + tracdata[0][1] * (dircos[0]);

              //cout << tracX << '\t' << tracY <<  '\t' << J << '\t' << Jmod << endl;
              //cout << endl;

              //cout << dircos[0] << '\t' << dircos[1] << endl;


              for(ii=0;ii<=p;ii++)
              {
                 index = 3*ii;
                 fact1 = NN[ii] * Jmod;
                 resi[index]   += tracX * fact1 ;
                 resi[index+1] += tracY * fact1 ;
              }
           }
        }

        // side #2
        if(!CompareDoubles(tracdata[1][0],7777.0) || !CompareDoubles(tracdata[1][1],7777.0))
        {
           double NN[q+1];

           for(gp=0;gp<nGP2;gp++)   // loop over Gauss points
           {
              if( finite && followerLoadFlag )
                 NurbsShapeFunctions2DAlg2(surf1, startindex[0], startindex[1], 1.0, gausspoints2[gp], NN, J, dircos);
              else
                 NurbsShapeFunctions2DAlg2(surf0, startindex[0], startindex[1], 1.0, gausspoints2[gp], NN, J, dircos);

              Jmod = J * gaussweights2[gp] * thick;

              //printf("\t J, Jmod  = %12.6f \t %12.6f \n",J, Jmod);

              //   for axisymmetric problems compute radius
              if(axsy)
              {
                 if(finite && followerLoadFlag )
                    rad  = surf1->SurfacePoint(uvalues[1], knotsAtGPs[2*(nGP1*(gp+1)-1)+1]).CalcEuclid().x;
                 else
                    rad  = surf0->SurfacePoint(uvalues[1], knotsAtGPs[2*(nGP1*(gp+1)-1)+1]).CalcEuclid().x;

                 Jmod *= (twoPI * rad);
              }

              tracX = tracdata[1][0] * (-dircos[0]) + tracdata[1][1] * (-dircos[1]);
              tracY = tracdata[1][0] * (-dircos[1]) + tracdata[1][1] * (dircos[0]);

              //printf(" tracX ... %12.6f \t %12.6f \n", tracdata[1][0], tracdata[1][1]);
              //printf(" tracX ... %12.6f \t %12.6f \n", tracX, tracY);

              for(ii=0;ii<=q;ii++)
              {
                 index = 3 * ((p + 1) * (ii+1) - 1);
                 fact1 = NN[ii] * Jmod;
                 resi[index]   += tracX * fact1 ;
                 resi[index+1] += tracY * fact1 ;
              }
           }
        }

        // side #3
        if(!CompareDoubles(tracdata[2][0],7777.0) || !CompareDoubles(tracdata[2][1],7777.0))
        {
           double NN[p+1];

           for(gp=0;gp<nGP1;gp++)   // loop over Gauss points
           {
              if( finite && followerLoadFlag )
                 NurbsShapeFunctions2DAlg2(surf1, startindex[0], startindex[1], gausspoints1[gp], 1.0, NN, J, dircos);
              else
                 NurbsShapeFunctions2DAlg2(surf0, startindex[0], startindex[1], gausspoints1[gp], 1.0, NN, J, dircos);

              Jmod = J * gaussweights1[gp] * thick;

              //   for axisymmetric problems compute radius
              if(axsy)
              {
                 if(finite && followerLoadFlag )
                    rad  = surf1->SurfacePoint(knotsAtGPs[2*((nGP2-1)*nGP1+gp)], vvalues[1]).CalcEuclid().x;
                 else
                    rad  = surf0->SurfacePoint(knotsAtGPs[2*((nGP2-1)*nGP1+gp)], vvalues[1]).CalcEuclid().x;

                 Jmod *= (twoPI * rad);
              }

              tracX = -tracdata[2][0] * (-dircos[0]) + -tracdata[2][1] * (-dircos[1]);
              tracY = -tracdata[2][0] * (-dircos[1]) + -tracdata[2][1] * (dircos[0]);

              //printf(" tracX ... %12.6f \t %12.6f \n", tracdata[2][0], tracdata[2][1]);
              //printf(" tracX ... %12.6f \t %12.6f \n", tracX, tracY);

              for(ii=0;ii<=p;ii++)
              {
                 index = 3 * ((p+1) * q + ii);
                 fact1 = NN[ii] * Jmod;
                 resi[index]   += tracX * fact1 ;
                 resi[index+1] += tracY * fact1 ;
              }

           }
        }

        // side #4
        if(!CompareDoubles(tracdata[3][0],7777.0) || !CompareDoubles(tracdata[3][1],7777.0))
        {
           double NN[q+1];

           for(gp=0;gp<nGP2;gp++)   // loop over Gauss points
           {
              if( finite && followerLoadFlag )
                 NurbsShapeFunctions2DAlg2(surf1, startindex[0], startindex[1], -1.0, gausspoints2[gp],  NN, J, dircos);
              else
                 NurbsShapeFunctions2DAlg2(surf0, startindex[0], startindex[1], -1.0, gausspoints2[gp],  NN, J, dircos);

              Jmod = J * gaussweights2[gp] * thick;

              //   for axisymmetric problems compute radius
              if(axsy)
              {
                 if(finite && followerLoadFlag )
                    rad  = surf1->SurfacePoint(uvalues[0], knotsAtGPs[gp*nGP1*2+1]).CalcEuclid().x;
                 else
                    rad  = surf0->SurfacePoint(uvalues[0], knotsAtGPs[gp*nGP1*2+1]).CalcEuclid().x;
                 
                 Jmod *= (twoPI * rad);
              }

              tracX = -tracdata[3][0] * (-dircos[0]) + -tracdata[3][1] * (-dircos[1]);
              tracY = -tracdata[3][0] * (-dircos[1]) + -tracdata[3][1] * (dircos[0]);
              
              //cout << tracX << '\t' << tracY << '\t' << rad << '\t' << J << '\t' << Jmod << endl;

              for(ii=0;ii<=q;ii++)
              {
                 index = 3 * ((p+1) * ii);
                 fact1 = NN[ii] * Jmod;
                 resi[index]   += tracX * fact1;
                 resi[index+1] += tracY * fact1;
              }

           }
        }
    }
//cout << " AAAAAAAAAA " << endl;
//printForceVector();
*/
  return 0;
}
//


void NurbsElem2DStructSolidLSFEM3dof::discreteContourplot(int vartype, int varindex, int index, int nCol, double umin, double umax)
{
  if(index > nivGP)
  {
     cout << '\t' << " Error in NurbsElem2DStructSolidLSFEM3dof::contourplot " << endl;
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

/*
  int tt=0;
  for(int jj=0;jj<nGP2;jj++)
  {
     for(int ii=0;ii<nGP1;ii++)
     {
        cout << '\t' << outval[tt];
        tt++;
     }
  }

  cout << endl;
  cout << endl;
*/

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

  double x1[2], x2[2], x3[2], x4[2], u1;

  int ind1, ind2, nGP1p1;

  nGP1p1 = nGP1 + 1;

  count=0;
  for(jj=0;jj<nGP2;jj++)
  {
      ind1 = nGP1p1*jj;
      ind2 = nGP1p1*(jj+1);

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











void NurbsElem2DStructSolidLSFEM3dof::projectToKnots(bool extrapolateFlag, int vartype, int varindex, int index)
{
/*
   vals2project[0] = intVar2[indx];
   vals2project[1] = intVar2[(nGP1-1)*nivGP+indx];
   vals2project[2] = intVar2[nGP1*(nGP2-1)*nivGP+indx];
   vals2project[3] = intVar2[(nGP1*nGP2-1)*nivGP+indx];
*/

   //double outval[nGP];
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




void NurbsElem2DStructSolidLSFEM3dof::projectStress(int varindex, double* outval)
{
/*
   if(varindex > 5)
   {
       cout << '\t' << "    NurbsElem2DStructSolidLSFEM3dof::projectStress .... : Error in 'varindex' " << endl;
       return;
   }

   double F[4], detF=0.0, F33=0.0, stre[4], cc[4][4], N[nlbf], dN_dx[nlbf], dN_dy[nlbf];

   double rad0=0.0, rad=0.0, Jac, dt, b1, b2, b3, pres, fact;

   int  err, isw, count, count1, ll, ii, jj, index, gp1, gp2;

   int nivEL = nGP * nivGP;
   for(ii=0;ii<nivEL;ii++)
     intVar2[ii] = intVar1[ii];

   count = 1;   ll = 0;   err = 0;   isw = 3;
   dt = mpapTime.dt;

if(!finite)
{
if(ndof == 2)
{
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

        if(varindex < 4)
           outval[count1] = stre[varindex];
        else if(varindex == 4)
           outval[count1] = sqrt((pow(stre[0]-stre[1],2.0) + pow(stre[1]-stre[2], 2.0) + pow(stre[2]-stre[0], 2.0) + 6* stre[3]*stre[3])/2);
        else if(varindex == 5)
           outval[count1] = pres;

        count++;
        count1++;
        ll += nivGP;

  }//gp1
  }//gp2
}
else
{
   double  *values1 = &(surf1->Values[0][0]);
   double  *values2 = &(surf1->Values[1][0]);
   double  *values3 = &(surf1->Values[2][0]);

   int *tt = &(surf0->IEN[elenum][0]);

   double  BULK, mu, lambda;
        
   BULK = matDat[0];
   mu = matDat[1];
   lambda = BULK - 2.0*mu/3.0;

   count1 = 0;
   for(gp2=0;gp2<nGP2;gp2++)
   {
   for(gp1=0;gp1<nGP1;gp1++)
   {
        index = count1*2;

        surf0->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, Jac);

        surf1->deformationGradient(startindex[0], startindex[1], 1, dN_dx, dN_dy, F, detF);

        F33 = 1.0;

        //matlib2d_(matDat, F, &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);

        pres = 0.0;
        for(ii=0;ii<nlbf;ii++)
        {
           index = tt[ii];

           b1 = values1[index];
           b2 = values2[index];
           b3 = values3[index];

           pres += b3*N[ii];
        }

          //pres = surf2->computeValue(1, knotsAtGPs[index], knotsAtGPs[index+1]);

          //stre[0] = 2.0*mu*( 2.0*(F[0]-1.0) - 1.0*(F[3]-1.0))/3.0 + pres;
          //stre[1] = 2.0*mu*(-1.0*(F[0]-1.0) + 2.0*(F[3]-1.0))/3.0 + pres;
          //stre[2] = 2.0*mu*(-1.0*(F[0]-1.0) - 1.0*(F[3]-1.0))/3.0 + pres;

          stre[0] = mu*(F[0] - F[3]) + pres;
          stre[1] = mu*(F[3] - F[0]) + pres;
          stre[2] = pres;
          stre[3] = mu*(F[1]+F[2]);


        if(varindex < 4)
           outval[count1] = stre[varindex];
        else if(varindex == 4)
           outval[count1] = sqrt((pow(stre[0]-stre[1],2.0) + pow(stre[1]-stre[2], 2.0) + pow(stre[2]-stre[0], 2.0) + 6* stre[3]*stre[3])/2);
        else if(varindex == 5)
           outval[count1] = pres;

        count++;
        count1++;
        ll += nivGP;

  }//gp1
  }//gp2
}
}
else
{
   double  *values1 = &(surf1->Values[0][0]);
   double  *values2 = &(surf1->Values[1][0]);
   double  *values3 = &(surf1->Values[2][0]);

   int *tt = &(surf0->IEN[elenum][0]);

   double  BULK, mu, aa, bb, fact, J;

   MatrixXd  b(3,3), FF(3,3);
        
   BULK = matDat[0];
   mu = matDat[1];

   aa = matDat[4];
   bb = matDat[5];

   count1 = 0;
   for(gp2=0;gp2<nGP2;gp2++)
   {
   for(gp1=0;gp1<nGP1;gp1++)
   {
        index = count1*2;

        surf0->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, Jac);

        surf1->deformationGradient2(tt, N, dN_dx, dN_dy, FF, b, pres);

        J = FF(0,0)*FF(1,1) - FF(0,1)*FF(1,0);

        if(CompareDoubles(bb, 3.0))
          b(2,2) = 1.0;
        if(CompareDoubles(bb, 2.0))
          b(2,2) = 0.0;

        fact = mu/pow(J,aa);

        b *= fact;

        fact = -b.trace()/bb + pres;

        stre[0] = b(0,0) + fact;
        stre[1] = b(1,1) + fact;
        stre[2] = b(2,2) + fact;
        stre[3] = b(0,1);


        if(varindex < 4)
           outval[count1] = stre[varindex];
        else if(varindex == 4)
           outval[count1] = sqrt((pow(stre[0]-stre[1],2.0) + pow(stre[1]-stre[2], 2.0) + pow(stre[2]-stre[0], 2.0) + 6* stre[3]*stre[3])/2);
        else if(varindex == 5)
           outval[count1] = pres;

        count++;
        count1++;
        ll += nivGP;

  }//gp1
  }//gp2
}
*/
  return;
}




void NurbsElem2DStructSolidLSFEM3dof::projectStrain(int vartype, int varindex, double* outval)
{
/*
  if(varindex > 5)
   {
       cout << '\t' << "    NurbsElem2DStructSolidLSFEM3dof::projectStress .... : Error in 'varindex' " << endl;
       return;
   }
   if(vartype > 3)
   {
       cout << '\t' << "    NurbsElem2DStructSolidLSFEM3dof::projectStress .... : Error in 'vartype' " << endl;
       return;
   }

   double F[4], detF=0.0, F33=0.0, stre[4], cc[4][4], N[nlbf], dN_dx[nlbf], dN_dy[nlbf];

   double rad0=0.0, rad=0.0, Jac, dt;

   int  err, isw, count, count1, ll, ii, jj, index, gp1, gp2;

   int nivEL = nGP * nivGP;
   for(ii=0;ii<nivEL;ii++)
     intVar2[ii] = intVar1[ii];

   bool intVarFlag = false;
   if( nivGP > 0)
      intVarFlag = true;

   count = 1;   ll = 0;   err = 0;   isw = 3;
   dt = mpapTime.dt;

   count1 = 0;
   for(gp2=0;gp2<nGP2;gp2++)
   {
   for(gp1=0;gp1<nGP1;gp1++)
   {
        index = count1*2;

        surf0->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, Jac);

        surf1->deformationGradient(startindex[0], startindex[1], 1, dN_dx, dN_dy, F, detF);

        if(detF < 0.0)
        {
           cerr << "        NurbsElem2DStructSolidLSFEM3dof::projectStress.......Negative DEFORMATION GRADIENT   " << endl;
           return ;
        }

        //   for axisymmetric problems compute radius
        if(axsy)
        {
           rad0 = surf0->SurfacePoint(knotsAtGPs[index], knotsAtGPs[index+1]).CalcEuclid().x;
           rad  = surf1->SurfacePoint(knotsAtGPs[index], knotsAtGPs[index+1]).CalcEuclid().x;
        }

        if(sss == 1)  // plane stress
        {
          if(finite)
            F33 = 1.0/sqrt(detF);
          else
            F33 = 3.0 - F[0] - F[3];
        }
        else if(sss == 2)    // plane strain
          F33 = 1.0;
        else //(sss == 3)    // axisymmetric
          F33 = rad/rad0;

        matlib2d_(matDat, F, &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, (Element*) this);

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

  }//gp1
  }//gp2
*/
  return;
}






void NurbsElem2DStructSolidLSFEM3dof::projectIntVar(int index, double* outval)
{
   int ind1, ii, jj;

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



int NurbsElem2DStructSolidLSFEM3dof::calcStiffnessMatrix(double dt)
{

  return 0;
}









int NurbsElem2DStructSolidLSFEM3dof::calcMassMatrix(int lumpInd, double dt)
{
   /*
   *  lumpInd = 1 --> consistent Mass matrix

   *          = 2 --> Row-sum Mass lumping
   *          = 3 --> proportional Mass lumping
   */

/*
   double fact1, fact2, fact, dvol0=0.0, rad0=0.0, J=0.0;

   VectorArray<double> NN;   NN.setDim(nlbf);


    mass_local.setDim(nsize);
    for(int ii=0;ii<nsize;ii++)
    {

       mass_local[ii].setDim(nsize);
       mass_local[ii].zero();

    }


       // loop over Gauss points

       for(int gp1=0;gp1<gausspoints1.n;gp1++)
       {
          for(int gp2=0;gp2<gausspoints2.n;gp2++)
          {

              NN.zero();

              NurbsShapeFunctions2DAlg3(surf0, startindex[0], startindex[1], gausspoints1[gp1], gausspoints2[gp2], NN, J);

              dvol0 = J * gaussweights[count1] * thick;

              //   for axisymmetric problems compute radius
              if(axsy)

              {
                 rad0 = 0.0 ;
                 int loc_num=0;
                 for(int jj=0; jj<=surf0->q; jj++)

                 {
                    for(int ii=0; ii<=surf0->p; ii++)
                    {
                       rad0 += (NN[loc_num] * surf0->Pw[startindex[0]+ii][startindex[1]+jj].CalcEuclid().x);

                       loc_num++;
                    }
                 }
                 dvol0 *= (twoPI * rad0);

              }

              for(int ii=0;ii<nlbf;ii++)
              {

                  fact = NN[ii] * dvol0;
                  for(int jj=0;jj<nlbf;jj++)
                  {
                     fact1 = fact * NN[jj];

                     mass_local[2*ii][2*jj]     += fact1;
                     mass_local[2*ii+1][2*jj+1] += fact1;
                  }
              }

          }
       }
*/
  return 0;
}





int NurbsElem2DStructSolidLSFEM3dof::calcOutput(double u1, double v1)
{
/*
   u1 = v1 = gausspoints[0];

   double F[4], detF=0.0, F33=0.0, fact, fact1, fact2, fact3, fact4,
          stre[4], cc[4][4], bc[2][4];

   double rad0=0.0, rad=0.0, dvol=0.0, dvol0=0.0, r1dF33=0.0, r1drad0=0.0, dt, J;

   int   err   = 0,
         isw   = 3,
         count = 1,
         ll    = 0,
         nivEl, ii, jj;

   nivEl = nivGP * nGP;

    for(ii=0;ii<nivEl;ii++)
       intVar2[ii] = intVar1[ii];

   double  N[nlbf], dN_dx[nlbf], dN_dy[nlbf];

        NurbsShapeFunctions2DAlg11(surf0, startindex[0], startindex[1], knotsAtGPs[index], knotsAtGPs[index+1], N, dN_dx, dN_dy, J);

        //   for axisymmetric problems compute radius
        if(axsy)
        {
           rad0 = surf0->SurfacePoint(u1, v1).CalcEuclid().x;
           rad  = surf1->SurfacePoint(u1, v1).CalcEuclid().x;

           r1drad0 = 1.0 / rad0;
        }

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
        else //(sss == 3)    // axisymmetric
          F33 = rad/rad0;

        dt = mpapTime.dt;
        matlib2d_(matDat, F, &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, (Element*) this);

        if(err !=0)
          return 1;
      //
         cout << "      stress vector ...:  " ;
         for(int ii=0;ii<4;ii++)
         {
           cout << fixed << '\t' << stre[ii];
         }
         cout << endl;
      //
*/
  return 0;
}


void NurbsElem2DStructSolidLSFEM3dof::toPostprocess(int vartype, int varindex, int type, SparseMatrixXd&  coeffMat, VectorXd& rhsVec)
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

   double  *values1 = &(surf1->Values[0][0]);
   double  *values2 = &(surf1->Values[1][0]);
   double  *values3 = &(surf1->Values[2][0]);

   int *tt = &(surf0->IEN[elenum][0]);

   double  BULK, mu, lambda, pres, b1, b2, b3, fact;
        
   BULK = matDat[0];
   mu = matDat[1];
   lambda = BULK - 2.0*mu/3.0;

   count1 = 0;
   for(gp2=0;gp2<nGP2;gp2++)
   {
   for(gp1=0;gp1<nGP1;gp1++)
   {
        index = count1*2;

        surf0->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), &(NN(0)), dN_dx, dN_dy, Jac);

        surf1->deformationGradient(startindex[0], startindex[1], 1, dN_dx, dN_dy, F, detF);

        F33 = 1.0;

        //matlib2d_(matDat, F, &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);

        pres = 0.0;
        for(ii=0;ii<nlbf;ii++)
        {
           index = tt[ii];

           b1 = values1[index];
           b2 = values2[index];
           b3 = values3[index];

           pres += b3*NN[ii];
        }

          //pres = surf2->computeValue(1, knotsAtGPs[index], knotsAtGPs[index+1]);

          fact = pres - (stre[0]+stre[1]+stre[2])/3.0 ;
          
          stre[0] += fact;
          stre[1] += fact;
          stre[2] += fact;

          //stre[0] = 2.0*mu*( 2.0*(F[0]-1.0) - 1.0*(F[3]-1.0))/3.0 + pres;
          //stre[1] = 2.0*mu*(-1.0*(F[0]-1.0) + 2.0*(F[3]-1.0))/3.0 + pres;
          //stre[2] = 2.0*mu*(-1.0*(F[0]-1.0) - 1.0*(F[3]-1.0))/3.0 + pres;

          stre[0] = mu*(F[0] - F[3]) + pres;
          stre[1] = mu*(F[3] - F[0]) + pres;
          stre[2] = pres;
          stre[3] = mu*(F[1]+F[2]);

        Nlocal += NN * NN.transpose();
        
        if(varindex < 4)
           rhslocal += NN * stre[varindex];
        else if(varindex == 4)
           rhslocal += NN * sqrt((pow(stre[0]-stre[1],2.0) + pow(stre[1]-stre[2], 2.0) + pow(stre[2]-stre[0], 2.0) + 6.0 * stre[3]*stre[3])/2.0);
        else if(varindex == 5)
           rhslocal += NN * pres;

        count++;
        count1++;
        ll += nivGP;

  }//gp1
  }//gp2

    for(ii=0;ii<nlbf;ii++)
    {
       row = tt[ii];
       rhsVec(row) += rhslocal[ii];
       for(jj=0;jj<nlbf;jj++)
       {
          col = tt[jj];
          coeffMat.coeffRef(row, col) += Nlocal(ii, jj);
       }
    }
*/

  return;
}







