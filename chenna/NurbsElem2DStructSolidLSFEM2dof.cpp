
#include <math.h>
#include "Debug.h"
#include "MpapTime.h"
#include "Plot.h"
#include "NurbsElem2DStructSolidLSFEM2dof.h"
#include "NurbsShapeFunctions.h"
#include <assert.h>
#include "ComputerTime.h"
#include "TimeFunction.h"


using namespace std;

extern ComputerTime       computerTime;
extern MpapTime mpapTime;
extern Plot plot;
extern List<TimeFunction> timeFunction;


NurbsElem2DStructSolidLSFEM2dof::NurbsElem2DStructSolidLSFEM2dof()
{
  if (debug) cout << " constructor NurbsElem2DStructSolidLSFEM2dof\n\n";
}



NurbsElem2DStructSolidLSFEM2dof::~NurbsElem2DStructSolidLSFEM2dof()
{
  if (debug) cout << " destructor NurbsElem2DStructSolidLSFEM2dof\n\n";
}



void NurbsElem2DStructSolidLSFEM2dof::createTractionDataVariable()
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



int NurbsElem2DStructSolidLSFEM2dof::calcStiffnessAndResidual()
{
  if(finite)
    calcStiffnessAndResidual2();
  else
    calcStiffnessAndResidual1();

   return 0;
}


int NurbsElem2DStructSolidLSFEM2dof::calcLoadVector()
{
  if(finite)
    calcLoadVector2();
  else
    calcLoadVector1();

   return 0;
}

int NurbsElem2DStructSolidLSFEM2dof::applyDirichletBCs()
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

     //printf(" \t %12.6f \t %12.6f \n", dispdata[ii][0], dispdata[ii][1]);
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

   //cout << " dispflag ... = " << dispflag << endl;

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





int NurbsElem2DStructSolidLSFEM2dof::calcStiffnessAndResidual1()
{
/*
   int ii, jj, gp1, gp2, TI, TIp1, count1, index;

   double  fact, dvol0, Jac, dt, ff1, ff2, lambda, BULK, mu, b1, b2, b3;

   double  N[nlbf], dN_dx[nlbf], d2N_dx2[nlbf], dN_dy[nlbf], d2N_dy2[nlbf], d2N_dxy[nlbf], d2N_dyx[nlbf];
   VectorXd  f(2);
   MatrixXd  D2u(2,3), D(nsize,2);
   
   BULK = matDat[0];
   mu = matDat[1];
   lambda = BULK - 2.0*mu/3.0;

   ff1 = 2.0*mu+lambda;
   ff2 = mu+lambda;
   
   bforce[0] = 0.0;
   bforce[1] = 1.0;
   rho0 = 1.0;

   double *gaussweights = &(surf0->gaussweights[0]);

   double  *values1 = &(surf1->Values[0][0]);
   double  *values2 = &(surf1->Values[1][0]);
   
   int *tt = &(surf0->IEN[elenum][0]);

   for(ii=0;ii<nsize;ii++)
     stiffness_local[ii].zero();
   
   resi.zero();

   count1 = 0;
   for(gp2=0;gp2<nGP2;gp2++)
   {
   for(gp1=0;gp1<nGP1;gp1++)
   {
        index = count1*2;

        surf0->ShapeFunDerivatives2(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, d2N_dx2, d2N_dy2, d2N_dxy, d2N_dyx, Jac);

        dvol0 = Jac * gaussweights[count1] * JacMultFact;
        count1++;

        //surf1->deformationGradient(startindex[0], startindex[1], 1, &dN_dx(0), &dN_dy(0), F, detF);

        //for(ii=0;ii<nlbf;ii++)
          //printf(" \t %14.8f\t %14.8f\t%14.8f\t %14.8f\t%14.8f\n", N[ii], dN_dx[ii], dN_dy[ii], d2N_dx2[ii], d2N_dy2[ii]);
        //printf("\n\n");

        f.setZero();
        D2u.setZero();

        f(0) = -bforce[0] * rho0 ;
        f(1) = -bforce[1] * rho0 ;

        for(ii=0;ii<nlbf;ii++)
        {
           TI = tt[ii];
             
           b1 = values1[TI];
           b2 = values2[TI];

           D2u(0,0) += b1*d2N_dx2[ii];
           D2u(0,1) += b1*d2N_dxy[ii];
           D2u(0,2) += b1*d2N_dy2[ii];

           D2u(1,0) += b2*d2N_dx2[ii];
           D2u(1,1) += b2*d2N_dxy[ii];
           D2u(1,2) += b2*d2N_dy2[ii];

           /////////////////////////////////////

           TI   =  2*ii;
           TIp1 =  TI+1;

           D(TI,0)   = ff1*d2N_dx2[ii] + mu*d2N_dy2[ii];
           D(TIp1,0) = ff2*d2N_dxy[ii];

           D(TI,1)   = D(TIp1,0);
           D(TIp1,1) = mu*d2N_dx2[ii] + ff1*d2N_dy2[ii];
        }

        //printf(" %14.8f \t%14.8f \t%14.8f \t%14.8f \t%14.8f \t%14.8f \t%14.8f \t%14.8f\n\n", ff, Du(0), Du(1), dp(0), dp(1), f(0), f(1), pres);

        f(0) -= (ff1*D2u(0,0) + mu*D2u(0,2) + ff2*D2u(1,1));
        f(1) -= (ff2*D2u(0,1) + mu*D2u(1,0) + ff1*D2u(1,2));

        for(ii=0;ii<nsize;ii++)
        {
           resi[ii] += dvol0*(D(ii,0)*f(0)+D(ii,1)*f(1));

           for(jj=0;jj<nsize;jj++)
             stiffness_local[ii][jj]  +=  dvol0*(D(ii,0)*D(jj,0)+D(ii,1)*D(jj,1)) ;
        }
    }//gp1
    }//gp2

//printStiffnessMatrix();
//printf("\n\n");
//printForceVector();
*/
   return 0;
}




int NurbsElem2DStructSolidLSFEM2dof::calcLoadVector1()
{
/*
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
   
      int *tt = &(surf0->IEN[elenum][0]);

      int p = surf0->p, q = surf0->q, ii, jj, gp, TI, TIp1, index, TJ, TJp1;
      double  NN[nlbf], dN_dx[nlbf], dN_dy[nlbf], Jac, fact, fact1, fact2, ALPHA, ALPHA1, ff1, ff2;
      double  J, Jmod, dircos[2], BULK, mu, b1, b2, b3, params[2], val1, val2, lambda;
        
      BULK = matDat[0];
      mu = matDat[1];
      lambda = BULK - 2.0*mu/3.0;

      ff1 = 2.0*mu+lambda;
      ff2 = mu+lambda;

      ALPHA = matDat[2];
      ALPHA1 = matDat[3];

      ALPHA1 = ALPHA = 1.0;

      bool pout = false;

      pout = (bool) matDat[3];

      // side #1
      if(!CompareDoubles(tracdata[0][0],7777) || !CompareDoubles(tracdata[0][1],7777))
      {
          double  N[p+1];

          val1 = 0.5*uvalues[2];
          val2 = 0.5*(uvalues[0]+uvalues[1]);

          params[1] = 0.0;
          for(gp=0;gp<nGP1;gp++)   // loop over Gauss points
          {
              params[0] = val1*gausspoints1[gp] + val2;
           
              NurbsShapeFunctions2DAlg2(surf0, startindex[0], startindex[1], gausspoints1[gp], -1.0, N, J, dircos);

              surf0->ShapeFunDerivatives(&(startindex[0]), params, NN, dN_dx, dN_dy, Jac);

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
              for(ii=0;ii<nlbf;ii++)
              {
                 index = tt[ii];

                 b1 = values1[index];
                 b2 = values2[index];

                 F(0,0) += b1 * dN_dx[ii];
                 F(0,1) += b1 * dN_dy[ii];
                 F(1,0) += b2 * dN_dx[ii];
                 F(1,1) += b2 * dN_dy[ii];

                 TI = 2*ii;
                 TIp1 = TI+1;

                 D(TI,0)   = ff1*dN_dx[ii]*dircos[0] + mu*dN_dy[ii]*dircos[1];
                 D(TIp1,0) = lambda*dN_dy[ii]*dircos[0] + mu*dN_dx[ii]*dircos[1];

                 D(TI,1)   = mu*dN_dy[ii]*dircos[0] + lambda*dN_dx[ii]*dircos[1];
                 D(TIp1,1) = mu*dN_dx[ii]*dircos[0] + ff1*dN_dy[ii]*dircos[1];
              }

              fact = mu*(F(0,1)+F(1,0));
              
              res(0) -= ((ff1*F(0,0) + lambda*F(1,1))*dircos[0] + fact*dircos[1]);
              res(1) -= (fact*dircos[0] + (lambda*F(0,0) + ff1*F(1,1))*dircos[1]);

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
           double   N[p+1];

           val1 = 0.5*uvalues[2];
           val2 = 0.5*(uvalues[0]+uvalues[1]);

           params[1] = 1.0;
           for(gp=0;gp<nGP1;gp++)   // loop over Gauss points
           {
              params[0] = val1*gausspoints1[gp] + val2;
           
              NurbsShapeFunctions2DAlg2(surf0, startindex[0], startindex[1], gausspoints1[gp], 1.0, N, J, dircos);

              surf0->ShapeFunDerivatives(&(startindex[0]), params, NN, dN_dx, dN_dy, Jac);

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
              for(ii=0;ii<nlbf;ii++)
              {
                 index = tt[ii];

                 b1 = values1[index];
                 b2 = values2[index];

                 F(0,0) += b1 * dN_dx[ii];
                 F(0,1) += b1 * dN_dy[ii];
                 F(1,0) += b2 * dN_dx[ii];
                 F(1,1) += b2 * dN_dy[ii];

                 TI = 2*ii;
                 TIp1 = TI+1;

                 D(TI,0)   = ff1*dN_dx[ii]*dircos[0] + mu*dN_dy[ii]*dircos[1];
                 D(TIp1,0) = lambda*dN_dy[ii]*dircos[0] + mu*dN_dx[ii]*dircos[1];

                 D(TI,1)   = mu*dN_dy[ii]*dircos[0] + lambda*dN_dx[ii]*dircos[1];
                 D(TIp1,1) = mu*dN_dx[ii]*dircos[0] + ff1*dN_dy[ii]*dircos[1];
              }

              fact = mu*(F(0,1)+F(1,0));
              
              res(0) -= ((ff1*F(0,0) + lambda*F(1,1))*dircos[0] + fact*dircos[1]);
              res(1) -= (fact*dircos[0] + (lambda*F(0,0) + ff1*F(1,1))*dircos[1]);

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
           double   N[q+1];

           val1 = 0.5*vvalues[2];
           val2 = 0.5*(vvalues[0]+vvalues[1]);

           params[0] = 1.0;
           for(gp=0;gp<nGP2;gp++)   // loop over Gauss points
           {
              params[1] = val1*gausspoints2[gp] + val2;
           
              NurbsShapeFunctions2DAlg2(surf0, startindex[0], startindex[1], 1.0, gausspoints2[gp], N, J, dircos);

              surf0->ShapeFunDerivatives(&(startindex[0]), params, NN, dN_dx, dN_dy, Jac);

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
              for(ii=0;ii<nlbf;ii++)
              {
                 index = tt[ii];

                 b1 = values1[index];
                 b2 = values2[index];

                 F(0,0) += b1 * dN_dx[ii];
                 F(0,1) += b1 * dN_dy[ii];
                 F(1,0) += b2 * dN_dx[ii];
                 F(1,1) += b2 * dN_dy[ii];

                 TI = 2*ii;
                 TIp1 = TI+1;

                 D(TI,0)   = ff1*dN_dx[ii]*dircos[0] + mu*dN_dy[ii]*dircos[1];
                 D(TIp1,0) = lambda*dN_dy[ii]*dircos[0] + mu*dN_dx[ii]*dircos[1];

                 D(TI,1)   = mu*dN_dy[ii]*dircos[0] + lambda*dN_dx[ii]*dircos[1];
                 D(TIp1,1) = mu*dN_dx[ii]*dircos[0] + ff1*dN_dy[ii]*dircos[1];
              }

              fact = mu*(F(0,1)+F(1,0));
              
              res(0) -= ((ff1*F(0,0) + lambda*F(1,1))*dircos[0] + fact*dircos[1]);
              res(1) -= (fact*dircos[0] + (lambda*F(0,0) + ff1*F(1,1))*dircos[1]);

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
           double   N[q+1];

           val1 = 0.5*vvalues[2];
           val2 = 0.5*(vvalues[0]+vvalues[1]);

           params[0] = 0.0;
           for(gp=0;gp<nGP2;gp++)   // loop over Gauss points
           {
              params[1] = val1*gausspoints2[gp] + val2;
           
              NurbsShapeFunctions2DAlg2(surf0, startindex[0], startindex[1], -1.0, gausspoints2[gp], N, J, dircos);

              surf0->ShapeFunDerivatives(&(startindex[0]), params, NN, dN_dx, dN_dy, Jac);

              Jmod = J * gaussweights2[gp] ;

              dircos[0] *= -1.0;
              dircos[1] *= -1.0;

              res(0) = tracdata[3][0] * (-dircos[0]) + tracdata[3][1] * (-dircos[1]);
              res(1) = tracdata[3][0] * (-dircos[1]) + tracdata[3][1] * (dircos[0]);

              //printf(" tracX and tracY ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n", res(0), res(1), J, Jmod, Jac);
              //printf(" dircos ... %12.6f \t %12.6f \n", dircos[0], dircos[1]);

              res.setZero();
              for(ii=0;ii<nlbf;ii++)
              {
                 index = tt[ii];

                 b1 = values1[index];
                 b2 = values2[index];

                 res(0) -= b1*NN[ii];
                 res(1) -= b2*NN[ii];
              }

              Jmod *= ALPHA;
              for(ii=0;ii<nlbf;ii++)
              {
                TI   = 2*ii;
                TIp1 = TI+1;

                fact1 = NN[ii]*Jmod;

                resi[TI]   += fact1*res(0);
                resi[TIp1] += fact1*res(1);

                for(jj=0;jj<nlbf;jj++)
                {
                  TJ   = 2*jj;
                  TJp1 = TJ+1;

                  fact = fact1 * NN[jj];

                  stiffness_local[TI][TJ]      +=  fact;
                  stiffness_local[TIp1][TJp1]  +=  fact;
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
*/
  return 0;
}




void NurbsElem2DStructSolidLSFEM2dof::discreteContourplot(int vartype, int varindex, int index, int nCol, double umin, double umax)
{
  if(index > nivGP)
  {
     cout << '\t' << " Error in NurbsElem2DStructSolidLSFEM2dof::contourplot " << endl;
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











void NurbsElem2DStructSolidLSFEM2dof::projectToKnots(bool extrapolateFlag, int vartype, int varindex, int index)
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




void NurbsElem2DStructSolidLSFEM2dof::projectStress(int varindex, double* outval)
{
/*
   if(varindex > 5)
   {
       cout << '\t' << "    NurbsElem2DStructSolidLSFEM2dof::projectStress .... : Error in 'varindex' " << endl;
       return;
   }

   double F[4], detF=0.0, F33=0.0, stre[4], cc[4][4], N[nlbf], dN_dx[nlbf], dN_dy[nlbf];
   double  Jac, dt, b1, b2, b3, pres, fact;

   int  err, isw, count, count1, ll, ii, jj, index, gp1, gp2;
   int nivEL = nGP * nivGP;

   for(ii=0;ii<nivEL;ii++)
     intVar2[ii] = intVar1[ii];

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
*/
  return;
}




void NurbsElem2DStructSolidLSFEM2dof::projectStrain(int vartype, int varindex, double* outval)
{
/*
   if(varindex > 5)
   {
       cout << '\t' << "    NurbsElem2DStructSolidLSFEM2dof::projectStress .... : Error in 'varindex' " << endl;
       return;
   }
   if(vartype > 3)
   {
       cout << '\t' << "    NurbsElem2DStructSolidLSFEM2dof::projectStress .... : Error in 'vartype' " << endl;
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
           cerr << "        NurbsElem2DStructSolidLSFEM2dof::projectStress.......Negative DEFORMATION GRADIENT   " << endl;
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






void NurbsElem2DStructSolidLSFEM2dof::projectIntVar(int index, double* outval)
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



int NurbsElem2DStructSolidLSFEM2dof::calcStiffnessMatrix(double dt)
{

  return 0;
}









int NurbsElem2DStructSolidLSFEM2dof::calcMassMatrix(int lumpInd, double dt)
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





int NurbsElem2DStructSolidLSFEM2dof::calcOutput(double u1, double v1)
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





void NurbsElem2DStructSolidLSFEM2dof::toPostprocess(int vartype, int varindex, int type, SparseMatrixXd&  coeffMat, VectorXd& rhsVec)
{
/*
   double F[4], detF=0.0, F33, Jac, dt, dN_dx[nlbf], dN_dy[nlbf], stre[4], cc[4][4];

   int   err,  isw,  count,  count1, index, ll = 0, ii, jj, gp1, gp2, row, col, *tt;

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

        surf0->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), &(NN(0)), dN_dx, dN_dy, Jac);

        surf1->deformationGradient(startindex[0], startindex[1], 1, dN_dx, dN_dy, F, detF);

        F33 = 1.0;    // plane strain

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



/*
void NurbsElem2DStructSolidLSFEM2dof::toPostprocess(int vartype, int varindex, int type, SparseMatrixXd&  coeffMat, VectorXd& rhsVec)
{

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

        surf0->ShapeFunDerivatives3(&(startindex[0]), &(knotsAtGPs[index]), &(NN(0)), dN_dx, dN_dy, Jac);

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


  return;
}
*/






