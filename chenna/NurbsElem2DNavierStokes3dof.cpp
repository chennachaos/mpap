
#include <math.h>
#include "Debug.h"
#include "MpapTime.h"
#include "Plot.h"
#include "NurbsElem2DNavierStokes3dof.h"
#include "NurbsShapeFunctions.h"
#include <assert.h>
#include "ComputerTime.h"
#include "Functions.h"

using namespace std;
using namespace Eigen;

extern ComputerTime       computerTime;
extern MpapTime mpapTime;
extern Plot plot;
extern List<TimeFunction> timeFunction;


NurbsElem2DNavierStokes3dof::NurbsElem2DNavierStokes3dof()
{
  if (debug) cout << " constructor NurbsElem2DNavierStokes3dof\n\n";
}



NurbsElem2DNavierStokes3dof::~NurbsElem2DNavierStokes3dof()
{
  if (debug) cout << " destructor NurbsElem2DNavierStokes3dof\n\n";
}




void NurbsElem2DNavierStokes3dof::createTractionDataVariable()
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



/*
int NurbsElem2DNavierStokes3dof::calcStiffnessAndResidual()
{
   double  fact, dvol0, Jac, nu, s, trgradu, b1, b2, b3, ci;
   double  N[nlbf], dN_dx[nlbf], d2N_dx2[nlbf], dN_dy[nlbf], d2N_dy2[nlbf];

   VectorXd  vel(2), f(3), dp(2), Du(2), res(3);
   MatrixXd  F(2,2), D(nsize,3);

   nu  =  1.0/elmDat[4];

   int   count1, ii, jj, gp1, gp2, index, TI, TIp1, TIp2, TJ, TJp1, TJp2;

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

        for(ii=0;ii<nlbf;ii++)
          printf(" \t %14.8f\t %14.8f\t%14.8f\t %14.8f\t%14.8f\n", N[ii], dN_dx[ii], dN_dy[ii], d2N_dx2[ii], d2N_dy2[ii]);
          
        f.setZero();
        res.setZero();
        vel.setZero();
        F.setZero();
        dp.setZero();
        Du.setZero();

        for(ii=0;ii<nlbf;ii++)
        {
           TI = tt[ii];
             
           b1 = values1[TI];
           b2 = values2[TI];
           b3 = values3[TI];
           
           vel(0) += b1 * N[ii];
           vel(1) += b2 * N[ii];
           
           fact = d2N_dx2[ii]+d2N_dy2[ii];

           Du(0)  += b1*fact;
           Du(1)  += b2*fact;
           
           dp(0)  += b3 * dN_dx[ii];
           dp(1)  += b3 * dN_dy[ii];

           F(0,0) += b1 * dN_dx[ii];
           F(0,1) += b1 * dN_dy[ii];
           F(1,0) += b2 * dN_dx[ii];
           F(1,1) += b2 * dN_dy[ii];
        }

        for(ii=0;ii<nlbf;ii++)
        {
           TI   =  3*ii;
           TIp1 =  TI+1;
           TIp2 =  TI+2;

           b1 = dN_dx[ii];
           b2 = dN_dy[ii];
           b3 = N[ii];

           ci = vel(0)*b1 + vel(1)*b2 - nu * (d2N_dx2[ii]+d2N_dy2[ii]);
             
           D(TI,0)    = ci + F(0,0)*b3;
           D(TIp1,0)  = F(0,1)*b3;
           D(TIp2,0)  = b1;

           D(TI,1)    = F(1,0)*b3;
           D(TIp1,1)  = ci + F(1,1)*b3;
           D(TIp2,1)  = b2;
           
           D(TI,2)    = b1;
           D(TIp1,2)  = b2;
           D(TIp2,2)  = 0.0;
        }

        res(0) = -F(0,0)*vel(0) - F(0,1)*vel(1) + nu*Du(0) - dp(0);
        res(1) = -F(1,0)*vel(0) - F(1,1)*vel(1) + nu*Du(1) - dp(1);
        res(2) = -F(0,0) - F(1,1);

        //f -= (F*vel - nu*Du + dp); //force due to linearisation of Navier-Stokes equations

        printf("\n\n \t %14.8f\t %14.8f\t%14.8f\t %14.8f\t%14.8f\t%14.8f\t%14.8f\t%14.8f\t%14.8f\n\n\n", f(0), f(1), F(0,0), F(0,1), F(1,0), F(1,1), JacMultFact, Jac, dvol0);

        for(ii=0;ii<nsize;ii++)
        {
           resi[ii] += dvol0*(D(ii,0)*res(0)+D(ii,1)*res(1)+D(ii,2)*res(2));

           for(jj=0;jj<nsize;jj++)
             stiffness_local[ii][jj]  +=  dvol0*(D(ii,0)*D(jj,0)+D(ii,1)*D(jj,1) + D(ii,2)*D(jj,2)) ;
        }

        //printf(" \t %14.8f \t%14.8f \t %14.8f\n", res(0), res(1), res(2));
  }//gp1
  }//gp2

//if(elenum == 0)
  //printStiffnessMatrix();

  //printForceVector();  
  //cout << endl;  cout << endl;

  return 0;
}
*/



int NurbsElem2DNavierStokes3dof::calcStiffnessAndResidual()
{
/*
   double  fact, dvol0, Jac, rho, mu, s, trgradu, b1, b2, b3, ci;
   double  N[nlbf], dN_dx[nlbf], d2N_dx2[nlbf], dN_dy[nlbf], d2N_dy2[nlbf], d2N_dx[nlbf], d2N_dyx[nlbf];

   VectorXd  vel(2), f(3), dp(2), Du(2), res(3);
   MatrixXd  F(2,2), D(nsize,3);

   rho = elmDat[3];
   mu = elmDat[4];

   int   count1, ii, jj, gp1, gp2, index, TI, TIp1, TIp2, TJ, TJp1, TJp2;

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

        surf0->ShapeFunDerivatives2(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, d2N_dx2, d2N_dy2, d2N_dx, d2N_dyx, Jac);

        dvol0 = Jac * gaussweights[count1] * JacMultFact;
        count1++;        

        //for(ii=0;ii<nlbf;ii++)
          //printf(" \t %14.8f\t %14.8f\t%14.8f\t %14.8f\t%14.8f\n", N[ii], dN_dx[ii], dN_dy[ii], d2N_dx2[ii], d2N_dy2[ii]);

        f.setZero();
        res.setZero();
        vel.setZero();
        F.setZero();
        dp.setZero();
        Du.setZero();

        for(ii=0;ii<nlbf;ii++)
        {
           TI = tt[ii];
             
           b1 = values1[TI];
           b2 = values2[TI];
           b3 = values3[TI];
           
           vel(0) += b1 * N[ii];
           vel(1) += b2 * N[ii];
           
           fact = d2N_dx2[ii]+d2N_dy2[ii];

           Du(0)  += b1*fact;
           Du(1)  += b2*fact;
           
           dp(0)  += b3 * dN_dx[ii];
           dp(1)  += b3 * dN_dy[ii];

           F(0,0) += b1 * dN_dx[ii];
           F(0,1) += b1 * dN_dy[ii];
           F(1,0) += b2 * dN_dx[ii];
           F(1,1) += b2 * dN_dy[ii];
        }

        for(ii=0;ii<nlbf;ii++)
        {
           TI   =  3*ii;
           TIp1 =  TI+1;
           TIp2 =  TI+2;

           b1 = dN_dx[ii];
           b2 = dN_dy[ii];
           b3 = N[ii];

           ci = rho*(vel(0)*b1 + vel(1)*b2) - mu*(d2N_dx2[ii]+d2N_dy2[ii]);

           D(TI,0)    = rho*F(0,0)*b3 + ci;
           D(TIp1,0)  = rho*F(0,1)*b3;
           D(TIp2,0)  = b1;

           D(TI,1)    = rho*F(1,0)*b3;
           D(TIp1,1)  = rho*F(1,1)*b3 + ci;
           D(TIp2,1)  = b2;
           
           D(TI,2)    = b1;
           D(TIp1,2)  = b2;
           D(TIp2,2)  = 0.0;
        }

        res(0) = rho*( F(0,0)*vel(0) + F(0,1)*vel(1) ) - mu*Du(0) + dp(0);
        res(1) = rho*( F(1,0)*vel(0) + F(1,1)*vel(1) ) - mu*Du(1) + dp(1);
        res(2) = F(0,0) + F(1,1);

        //f -= (F*vel - nu*Du + dp); //force due to linearisation of Navier-Stokes equations

        //printf("\n\n \t %14.8f\t %14.8f\t%14.8f\t %14.8f\t%14.8f\t%14.8f\t%14.8f\t%14.8f\t%14.8f\n\n\n", f(0), f(1), F(0,0), F(0,1), F(1,0), F(1,1), JacMultFact, Jac, dvol0);

        for(ii=0;ii<nsize;ii++)
        {
           resi[ii] -= dvol0*(D(ii,0)*res(0)+D(ii,1)*res(1)+D(ii,2)*res(2));

           for(jj=0;jj<nsize;jj++)
             stiffness_local[ii][jj]  +=  dvol0*(D(ii,0)*D(jj,0)+D(ii,1)*D(jj,1) + D(ii,2)*D(jj,2)) ;
        }

        //printf(" \t %14.8f \t%14.8f \t %14.8f\n", res(0), res(1), res(2));
  }//gp1
  }//gp2

//if(elenum == 0)
  //printStiffnessMatrix();

  //printForceVector();  
  //cout << endl;  cout << endl;
*/
  return 0;
}
//





int NurbsElem2DNavierStokes3dof::calcInternalForces()
{
/*
   double  fact, dvol0, Jac, nu, s, b1, b2, b3, ci, cj;

   double  N[nlbf], dN_dx[nlbf], d2N_dx2[nlbf], dN_dy[nlbf], d2N_dy2[nlbf], d2N_dx[nlbf], d2N_dyx[nlbf];

   VectorXd  vel(2), f(3), dp(3), Du(3);
   
   MatrixXd  F(2,2), D(nsize,3);

   nu  =  1.0/elmDat[4];

   int   count1, ii, jj, gp1, gp2, index, TI, TIp1, TIp2;

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

        surf0->ShapeFunDerivatives2(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, d2N_dx2, d2N_dy2, d2N_dx, d2N_dyx, Jac);

        fact = gaussweights[count1] * JacMultFact;
        count1++;
        dvol0 = Jac * fact;
        
        //for(ii=0;ii<nlbf;ii++)
          // printf(" \t %14.8f\t %14.8f\t%14.8f\t %14.8f\t%14.8f\n", N[ii], dN_dx[ii], dN_dy[ii], d2N_dx2[ii], d2N_dy2[ii]);
          
        f.setZero();
        vel.setZero();
        F.setZero();
        dp.setZero();
        Du.setZero();

        for(ii=0;ii<nlbf;ii++)
        {
           TI = tt[ii];
             
           b1 = values1[TI];
           b2 = values2[TI];
           b3 = values3[TI];
           
           vel(0) += b1 * N[ii];
           vel(1) += b2 * N[ii];
           
           fact = d2N_dx2[ii]+d2N_dy2[ii];

           Du(0)  += b1*fact;
           Du(1)  += b2*fact;
           
           dp(0)  += b3 * dN_dx[ii];
           dp(1)  += b3 * dN_dy[ii];

           F(0,0) += b1 * dN_dx[ii];
           F(0,1) += b1 * dN_dy[ii];
           F(1,0) += b2 * dN_dx[ii];
           F(1,1) += b2 * dN_dy[ii];
        }

        for(ii=0;ii<nlbf;ii++)
        {
           TI   =  3*ii;
           TIp1 =  TI+1;
           TIp2 =  TI+2;

           b1 = dN_dx[ii];
           b2 = dN_dy[ii];
           b3 = N[ii];

           ci = vel(0)*b1 + vel(1)*b2 - nu * (d2N_dx2[ii]+d2N_dy2[ii]);
             
           D(TI,0)    = ci + (1.0*F(0,0)+0.0*F(1,1))*b3;
           D(TIp1,0)  = F(0,1)*b3;
           D(TIp2,0)  = b1;

           D(TI,1)    = F(1,0)*b3;
           D(TIp1,1)  = ci + (0.0*F(0,0)+1.0*F(1,1))*b3;
           D(TIp2,1)  = b2;
           
           D(TI,2)    = b1;
           D(TIp1,2)  = b2;
           D(TIp2,2)  = 0.0;
        }

        f(0) -= (F(0,0)*vel(0) + F(0,1)*vel(1) - nu*Du(0) + dp(0));
        f(1) -= (F(1,0)*vel(0) + F(1,1)*vel(1) - nu*Du(1) + dp(1));
        f(2) = -(F(0,0) + F(1,1));

        for(ii=0;ii<nsize;ii++)
           resi[ii] += dvol0*(D(ii,0)*f(0)+D(ii,1)*f(1)+D(ii,2)*f(2));
  }//gp1
  }//gp2

  //printForceVector();  
  //cout << endl;  cout << endl;
   if(tracflag)
   {
      double *gausspoints1 = &(surf0->gausspoints1[0]);
      double *gausspoints2 = &(surf0->gausspoints2[0]);
      double *gaussweights1 = &(surf0->gaussweights1[0]);
      double *gaussweights2 = &(surf0->gaussweights2[0]);

      resi.zero();
    
      VectorXd  trac(2), res(2);
      MatrixXd  D(nsize, 2), F(2,2);
    
      double  *values1 = &(surf1->Values[0][0]);
      double  *values2 = &(surf1->Values[1][0]);
      double  *values3 = &(surf1->Values[2][0]);
   
      int *tt = &(surf0->IEN[elenum][0]);

      Kovasznay  kova1(elmDat[4]);

      kova1.SetPressure(1.310741966654558776639305506250821053981781005859375);

        EPOINT  EP;

        //cout << "       elem... : " << elenum << endl;
        //
        cout << tracdata[0][0] << '\t' << tracdata[0][1] << endl;
        cout << tracdata[1][0] << '\t' << tracdata[1][1] << endl;
        cout << tracdata[2][0] << '\t' << tracdata[2][1] << endl;
        cout << tracdata[3][0] << '\t' << tracdata[3][1] << endl;
        //

        double  J, Jmod, dircos[2], BULK, mu, a1, a2, b1, b2, b3, params[2], val1, val2, d, pres;
        
        int p = surf0->p, q = surf0->q, ii, jj, gp, TI, TIp1, TIp2, index, nlbf2, TJ, TJp1, TJp2;
        double   NN[nlbf], dN_dx[nlbf], dN_dy[nlbf], Jac, fact, fact1, fact2, ALPHA, ALPHA1;

        ALPHA = matDat[2];
        ALPHA1 = matDat[3];

        if( matDat[3] > 1.0)
          ALPHA /= (uvalues[2]*vvalues[2]);

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

              surf0->ShapeFunDerivatives(&(startindex[0]), params, NN, dN_dx, dN_dy, Jac);

              Jmod = J * gaussweights1[gp];

              EP = surf0->SurfacePoint(params[0], params[1]).CalcEuclid();

              res(0) = kova1.computeValue(0, EP.x, EP.y);
              res(1) = kova1.computeValue(1, EP.x, EP.y);

              //printf(" tracX and tracY ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n", res(0), res(1), J, Jmod, Jac);
              //printf(" dircos ... %12.6f \t %12.6f \n", dircos[0], dircos[1]);

              //pres = 0.0;
              for(ii=0;ii<nlbf;ii++)
              {
                 index = tt[ii];

                 b1 = values1[index];
                 b2 = values2[index];
                 //b3 = values3[index];

                 res(0) -= b1*NN[ii];
                 res(1) -= b2*NN[ii];
                 //pres -= b3*NN[ii];
              }

              Jmod *= ALPHA;

              for(ii=0;ii<nlbf;ii++)
              {
                TI   = ndof*ii;
                TIp1 = TI+1;

                fact1 = NN[ii]*Jmod;

                resi[TI]   += fact1*res(0);
                resi[TIp1] += fact1*res(1);
              }
          }
          //cout << " side1 done " << endl;
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

              surf0->ShapeFunDerivatives(&(startindex[0]), params, NN, dN_dx, dN_dy, Jac);

              Jmod = J * gaussweights1[gp] ;

              EP = surf0->SurfacePoint(params[0], params[1]).CalcEuclid();

              res(0) = kova1.computeValue(0, EP.x, EP.y);
              res(1) = kova1.computeValue(1, EP.x, EP.y);

              //printf(" tracX and tracY ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n", res(0), res(1), J, Jmod, Jac);
              //printf(" dircos ... %12.6f \t %12.6f \n", dircos[0], dircos[1]);

              //pres = 0.0;
              for(ii=0;ii<nlbf;ii++)
              {
                 index = tt[ii];

                 b1 = values1[index];
                 b2 = values2[index];
                 //b3 = values3[index];

                 res(0) -= b1*NN[ii];
                 res(1) -= b2*NN[ii];
                 //pres -= b3*NN[ii];
              }

              Jmod *= ALPHA;

              for(ii=0;ii<nlbf;ii++)
              {
                TI   = ndof*ii;
                TIp1 = TI+1;
                //TIp2 = TI+2;

                fact1 = NN[ii]*Jmod;

                resi[TI]   += fact1*res(0);
                resi[TIp1] += fact1*res(1);
              }
           }
           //cout << " side3 done " << endl;
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

              surf0->ShapeFunDerivatives(&(startindex[0]), params, NN, dN_dx, dN_dy, Jac);

              Jmod = J * gaussweights2[gp] ;
              
              EP = surf0->SurfacePoint(params[0], params[1]).CalcEuclid();

              res(0) = kova1.computeValue(0, EP.x, EP.y);
              res(1) = kova1.computeValue(1, EP.x, EP.y);

              //printf(" tracX and tracY ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n", res(0), res(1), J, Jmod, Jac);
              //printf(" dircos ... %12.6f \t %12.6f \n", dircos[0], dircos[1]);

              //pres = 0.0;
              for(ii=0;ii<nlbf;ii++)
              {
                 index = tt[ii];

                 b1 = values1[index];
                 b2 = values2[index];
                 //b3 = values3[index];

                 res(0) -= b1*NN[ii];
                 res(1) -= b2*NN[ii];
                 //pres -= b3*NN[ii];
              }

              Jmod *= ALPHA;

              for(ii=0;ii<nlbf;ii++)
              {
                TI   = ndof*ii;
                TIp1 = TI+1;
                //TIp2 = TI+2;

                fact1 = NN[ii]*Jmod;

                resi[TI]   += fact1*res(0);
                resi[TIp1] += fact1*res(1);
              }
           }
           //cout << " side2 done " << endl;
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

              surf0->ShapeFunDerivatives(&(startindex[0]), params, NN, dN_dx, dN_dy, Jac);

              Jmod = J * gaussweights2[gp] ;

              EP = surf0->SurfacePoint(params[0], params[1]).CalcEuclid();

              res(0) = kova1.computeValue(0, EP.x, EP.y);
              res(1) = kova1.computeValue(1, EP.x, EP.y);

              //printf(" tracX and tracY ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n", res(0), res(1), J, Jmod, Jac);
              //printf(" dircos ... %12.6f \t %12.6f \n", dircos[0], dircos[1]);

              //pres = 0.0;
              for(ii=0;ii<nlbf;ii++)
              {
                 index = tt[ii];

                 b1 = values1[index];
                 b2 = values2[index];
                 //b3 = values3[index];

                 res(0) -= b1*NN[ii];
                 res(1) -= b2*NN[ii];
                 //pres -= b3*NN[ii];
              }

              Jmod *= ALPHA;

              for(ii=0;ii<nlbf;ii++)
              {
                TI   = ndof*ii;
                TIp1 = TI+1;
                //TIp2 = TI+2;

                fact1 = NN[ii]*Jmod;

                resi[TI]   += fact1*res(0);
                resi[TIp1] += fact1*res(1);
              }
           }
           //cout << " side4 done " << endl;
        }
//printForceVector();
    if(elenum == 0)
    {
       nlbf2 = q+1;
       double   N[nlbf2];

       params[0] = 0.0;
       params[1] = 0.0;
           
       NurbsShapeFunctions2DAlg2(surf0, startindex[0], startindex[1], -1.0, 0.0, N, J, dircos);

       surf0->ShapeFunDerivatives(&(startindex[0]), params, NN, dN_dx, dN_dy, Jac);

       Jmod = J ;

       //EP = surf0->SurfacePoint(params[0], params[1]).CalcEuclid();
       //kova1.computeVelocity(EP.x, EP.y, &(res(0)));

       //printf(" tracX and tracY ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n", res(0), res(1), J, Jmod, Jac);
       //printf(" dircos ... %12.6f \t %12.6f \n", dircos[0], dircos[1]);

       //pres = 0.0;

       EP = surf0->SurfacePoint(params[0], params[1]).CalcEuclid();

       pres = kova1.computeValue(2, EP.x, EP.y);

       pres -= values3[0] ;

       //printf("\t EP.x = %12.6f \t EP.y = %12.6f \n", EP.x, EP.y);
       //printf("\t Jmod = %12.6f \t pres = %12.6E \n", Jmod, pres);

       Jmod *= ALPHA;

       //Jmod *= 1000.0;

       resi[2] += pres * Jmod;
    }
    }
*/
  return 0;
}




int NurbsElem2DNavierStokes3dof::calcError(int index)
{
/*
    Kovasznay  kova1(elmDat[4]);

    //kova1.SetPressure(1.310741966654558776639305506250821053981781005859375);
    kova1.SetPressure(0.0);

    EPOINT  EP;

   int ii, jj, gp1, gp2, count1,  TI;

   double  *gaussweights = &(surf0->gaussweights[0]);
   double  *values1 = &(surf1->Values[0][0]);
   double  *values2 = &(surf1->Values[1][0]);
   double  *values3 = &(surf1->Values[2][0]);
   
   int *tt = &(surf0->IEN[elenum][0]);
   
   double  Jac, dvol0, fact, b1, b2, pres;

   double  N[nlbf], dN_dx[nlbf], dN_dy[nlbf], v[2];
    
   elemError = 0.0;

   if(index == 1) // x - velocity
   {
      count1 = 0;
      for(gp2=0;gp2<nGP2;gp2++)
      {
      for(gp1=0;gp1<nGP1;gp1++)
      {
          index = count1*2;

          surf0->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, Jac);

          EP = surf0->SurfacePoint(knotsAtGPs[index], knotsAtGPs[index+1]).CalcEuclid();

          dvol0 = Jac * gaussweights[count1] * JacMultFact;
          count1++;

          v[0] = kova1.computeValue(0, EP.x, EP.y);
          v[1] = kova1.computeValue(1, EP.x, EP.y);

          for(ii=0;ii<nlbf;ii++)
             v[0] -= values1[tt[ii]]* N[ii];

          fact = v[0]*v[0];
          
          elemError += ( fact * dvol0 );
      }
      }
   }
   if(index == 2) // y - velocity
   {
      count1 = 0;
      for(gp2=0;gp2<nGP2;gp2++)
      {
      for(gp1=0;gp1<nGP1;gp1++)
      {
          index = count1*2;
          surf0->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, Jac);

          EP = surf0->SurfacePoint(knotsAtGPs[index], knotsAtGPs[index+1]).CalcEuclid();

          dvol0 = Jac * gaussweights[count1] * JacMultFact;
          count1++;

          v[0] = kova1.computeValue(0, EP.x, EP.y);
          v[1] = kova1.computeValue(1, EP.x, EP.y);

          for(ii=0;ii<nlbf;ii++)
            v[1] -= values2[tt[ii]]* N[ii];

          fact = v[1]*v[1];
          
          elemError += ( fact * dvol0 );
      }
      }
      //printf("\n\n"); 
   }
   if(index == 3) // pressure
   {
      count1 = 0;
      for(gp2=0;gp2<nGP2;gp2++)
      {
      for(gp1=0;gp1<nGP1;gp1++)
      {
          index = count1*2;

          surf0->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, Jac);

          EP = surf0->SurfacePoint(knotsAtGPs[index], knotsAtGPs[index+1]).CalcEuclid();

          dvol0 = Jac * gaussweights[count1] * JacMultFact;
          count1++;

          pres = kova1.computeValue(2, EP.x, EP.y);

          for(ii=0;ii<nlbf;ii++)
            pres -= values3[tt[ii]] * N[ii];

          //printf("\t pres = %12.6E \n", pres);

          fact = pres*pres;
          
          elemError += ( fact * dvol0 );
      }
      }
      //printf("\n\n");
   }
   if(index == 4) // H1 norm in velocity
   {
      MatrixXd  F(2,2);

      count1 = 0;
      for(gp2=0;gp2<nGP2;gp2++)
      {
      for(gp1=0;gp1<nGP1;gp1++)
      {
          index = count1*2;

          surf0->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, Jac);

          EP = surf0->SurfacePoint(knotsAtGPs[index], knotsAtGPs[index+1]).CalcEuclid();

          dvol0 = Jac * gaussweights[count1] * JacMultFact;
          count1++;

          F.setZero();
          v[0] = kova1.computeValue(0, EP.x, EP.y);
          v[1] = kova1.computeValue(1, EP.x, EP.y);
          kova1.computeDerivatives(EP.x, EP.y, &(F(0,0)));

          for(ii=0;ii<nlbf;ii++)
          {
             TI = tt[ii];
             
             b1 = values1[TI];
             b2 = values2[TI];
           
             v[0] -= b1 * N[ii];
             v[1] -= b2 * N[ii];

             F(0,0) -= b1 * dN_dx[ii];
             F(0,1) -= b1 * dN_dy[ii];
             F(1,0) -= b2 * dN_dx[ii];
             F(1,1) -= b2 * dN_dy[ii];
          }

          fact = v[0]*v[0]+v[1]*v[1] + F(0,0)*F(0,0)+F(0,1)*F(0,1)+F(1,0)*F(1,0)+F(1,1)*F(1,1);
          
          //printf(" difference \t %12.8f \t %12.8f \n", v2(0), v2(1));
          
          elemError += ( fact * dvol0 );
      }
      }
   }
   if(index == 5) // LS functional
   {
       double  d2N_dx2[nlbf], d2N_dy2[nlbf], d2N_dxy[nlbf], d2N_dyx[nlbf];
       VectorXd  vel(2), f(3), dp(2), Du(2);
   
       MatrixXd  F(2,2);

       double  nu  =  1.0/elmDat[4], b3;

       int   TI, TIp1, TIp2;

       count1 = 0;
       for(gp2=0;gp2<nGP2;gp2++)
       {
       for(gp1=0;gp1<nGP1;gp1++)
       {
          index = count1*2;

          surf0->ShapeFunDerivatives2(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, d2N_dx2, d2N_dy2, d2N_dxy, d2N_dyx, Jac);

          dvol0 = Jac * gaussweights[count1] * JacMultFact;
          count1++;
        
          //for(ii=0;ii<nlbf;ii++)
            // printf(" \t %14.8f\t %14.8f\t%14.8f\t %14.8f\t%14.8f\n", N[ii], dN_dx[ii], dN_dy[ii], d2N_dx2[ii], d2N_dy2[ii]);
          
          f.setZero();
          vel.setZero();
          F.setZero();
          dp.setZero();
          Du.setZero();

          for(ii=0;ii<nlbf;ii++)
          {
             TI = tt[ii];
             
             b1 = values1[TI];
             b2 = values2[TI];
             b3 = values3[TI];
           
             vel(0) += b1 * N[ii];
             vel(1) += b2 * N[ii];
           
             fact = d2N_dx2[ii]+d2N_dy2[ii];

             Du(0)  += b1*fact;
             Du(1)  += b2*fact;
           
             dp(0)  += b3 * dN_dx[ii];
             dp(1)  += b3 * dN_dy[ii];

             F(0,0) += b1 * dN_dx[ii];
             F(0,1) += b1 * dN_dy[ii];
             F(1,0) += b2 * dN_dx[ii];
             F(1,1) += b2 * dN_dy[ii];
          }

          f(0) -= (F(0,0)*vel(0) + F(0,1)*vel(1) - nu*Du(0) + dp(0));
          f(1) -= (F(1,0)*vel(0) + F(1,1)*vel(1) - nu*Du(1) + dp(1));
          f(2) = -(F(0,0) + F(1,1));

          elemError += (f.norm()*dvol0);

      }//gp1
      }//gp2
   }
      //elemError = sqrt(elemError);///volume;

   if( (int) matDat[4] == 1)
     printf(" \t element = %5d ... \t ... elemError  =   %12.6E \n " , elenum, elemError);
*/
   return 0;
}




int NurbsElem2DNavierStokes3dof::calcLoadVector()
{
/*
   if(tracflag)
   {
      double *gausspoints1 = &(surf0->gausspoints1[0]);
      double *gausspoints2 = &(surf0->gausspoints2[0]);
      double *gaussweights1 = &(surf0->gaussweights1[0]);
      double *gaussweights2 = &(surf0->gaussweights2[0]);

      for(int ii=0;ii<nsize;ii++)
        stiffness_local[ii].zero();

      resi.zero();
    
      VectorXd  trac(2), res(2);
      MatrixXd  D(nsize, 2), F(2,2);
    
      double  *values1 = &(surf1->Values[0][0]);
      double  *values2 = &(surf1->Values[1][0]);
      double  *values3 = &(surf1->Values[2][0]);
   
      int *tt = &(surf0->IEN[elenum][0]);

      Kovasznay  kova1(elmDat[4]);

      //kova1.SetPressure(1.310741966654558776639305506250821053981781005859375);
      kova1.SetPressure(0.0);

        EPOINT  EP;

        //cout << "       elem... : " << elenum << endl;
        //
        cout << tracdata[0][0] << '\t' << tracdata[0][1] << endl;
        cout << tracdata[1][0] << '\t' << tracdata[1][1] << endl;
        cout << tracdata[2][0] << '\t' << tracdata[2][1] << endl;
        cout << tracdata[3][0] << '\t' << tracdata[3][1] << endl;
        //

        double  J, Jmod, dircos[2], BULK, mu, a1, a2, b1, b2, b3, params[2], val1, val2, d, pres;
        
        int p = surf0->p, q = surf0->q, ii, jj, gp, TI, TIp1, TIp2, index, nlbf2, TJ, TJp1, TJp2;
        double   NN[nlbf], dN_dx[nlbf], dN_dy[nlbf], Jac, fact, fact1, fact2, ALPHA, ALPHA1;

        ALPHA = matDat[3];
        ALPHA1 = matDat[4];

        if( matDat[3] > 1.0)
          ALPHA /= (uvalues[2]*vvalues[2]);

        //ALPHA = ALPHA1  = 1000;
        ALPHA1 = ALPHA = 1.0;
        //ALPHA1 = ALPHA;

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

              surf0->ShapeFunDerivatives(&(startindex[0]), params, NN, dN_dx, dN_dy, Jac);

              Jmod = J * gaussweights1[gp];

              EP = surf0->SurfacePoint(params[0], params[1]).CalcEuclid();

              res(0) = tracdata[0][0] *timeFunction[0].prop;
              res(1) = tracdata[0][1] *timeFunction[0].prop;

              //res[0] = kova1.computeValue(0, EP.x, EP.y);
              //res[1] = kova1.computeValue(1, EP.x, EP.y);

              //for(ii=0;ii<nlbf;ii++)
                //printf(" \t %14.8f \n", NN[ii]);

              //printf("\n\n tracX and tracY ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n\n\n", EP.x, EP.y, res(0), res(1), J, Jmod, Jac);
              //printf(" dircos ... %12.6f \t %12.6f \n", dircos[0], dircos[1]);

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
                TI   = ndof*ii;
                TIp1 = TI+1;

                fact1 = NN[ii]*Jmod;

                resi[TI]   += fact1*res(0);
                resi[TIp1] += fact1*res(1);

                for(jj=0;jj<nlbf;jj++)
                {
                  TJ   = ndof*jj;
                  TJp1 = TJ+1;

                  fact = fact1 * NN[jj];

                  stiffness_local[TI][TJ]      +=  fact;
                  stiffness_local[TIp1][TJp1]  +=  fact;
                }
              }
          }
          //cout << " side1 done " << endl;
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

              surf0->ShapeFunDerivatives(&(startindex[0]), params, NN, dN_dx, dN_dy, Jac);

              Jmod = J * gaussweights1[gp] ;

              EP = surf0->SurfacePoint(params[0], params[1]).CalcEuclid();

              res(0) = tracdata[2][0] *timeFunction[0].prop;
              res(1) = tracdata[2][1] *timeFunction[0].prop;

              //res[0] = kova1.computeValue(0, EP.x, EP.y);
              //res[1] = kova1.computeValue(1, EP.x, EP.y);

              //for(ii=0;ii<nlbf;ii++)
                //printf(" \t %14.8f \n", NN[ii]);

              //printf("\n\n tracX and tracY ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n\n\n", EP.x, EP.y, res(0), res(1), J, Jmod, Jac);


              //printf(" tracX and tracY ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n", res(0), res(1), J, Jmod, Jac);
              //printf(" dircos ... %12.6f \t %12.6f \n", dircos[0], dircos[1]);

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
                TI   = ndof*ii;
                TIp1 = TI+1;

                fact1 = NN[ii]*Jmod;

                resi[TI]   += fact1*res(0);
                resi[TIp1] += fact1*res(1);

                for(jj=0;jj<nlbf;jj++)
                {
                  TJ   = ndof*jj;
                  TJp1 = TJ+1;

                  fact = fact1 * NN[jj];

                  stiffness_local[TI][TJ]      +=  fact;
                  stiffness_local[TIp1][TJp1]  +=  fact;
                }
              }
           }
           //cout << " side3 done " << endl;
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

              surf0->ShapeFunDerivatives(&(startindex[0]), params, NN, dN_dx, dN_dy, Jac);

              Jmod = J * gaussweights2[gp] ;
              
              EP = surf0->SurfacePoint(params[0], params[1]).CalcEuclid();

              res(0) = tracdata[1][0] *timeFunction[0].prop;
              res(1) = tracdata[1][1] *timeFunction[0].prop;

              //res[0] = kova1.computeValue(0, EP.x, EP.y);
              //res[1] = kova1.computeValue(1, EP.x, EP.y);

              //for(ii=0;ii<nlbf;ii++)
                //printf(" \t %14.8f \n", NN[ii]);

              //printf("\n\n tracX and tracY ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n\n\n", EP.x, EP.y, res(0), res(1), J, Jmod, Jac);

              //printf(" tracX and tracY ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n", res(0), res(1), J, Jmod, Jac);
              //printf(" dircos ... %12.6f \t %12.6f \n", dircos[0], dircos[1]);

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
                TI   = ndof*ii;
                TIp1 = TI+1;

                fact1 = NN[ii]*Jmod;

                resi[TI]   += fact1*res(0);
                resi[TIp1] += fact1*res(1);

                for(jj=0;jj<nlbf;jj++)
                {
                  TJ   = ndof*jj;
                  TJp1 = TJ+1;

                  fact = fact1 * NN[jj];

                  stiffness_local[TI][TJ]      +=  fact;
                  stiffness_local[TIp1][TJp1]  +=  fact;
                }
              }
           }
           //cout << " side2 done " << endl;
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

              surf0->ShapeFunDerivatives(&(startindex[0]), params, NN, dN_dx, dN_dy, Jac);

              Jmod = J * gaussweights2[gp] ;

              EP = surf0->SurfacePoint(params[0], params[1]).CalcEuclid();

              res(0) = tracdata[3][0] *timeFunction[0].prop;
              res(1) = tracdata[3][1] *timeFunction[0].prop;

              //res[0] = kova1.computeValue(0, EP.x, EP.y);
              //res[1] = kova1.computeValue(1, EP.x, EP.y);

              //for(ii=0;ii<nlbf;ii++)
                //printf(" \t %14.8f \n", NN[ii]);

              //printf("\n\n tracX and tracY ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n\n\n", EP.x, EP.y, res(0), res(1), J, Jmod, Jac);


              //printf(" tracX and tracY ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n", res(0), res(1), J, Jmod, Jac);
              //printf(" dircos ... %12.6f \t %12.6f \n", dircos[0], dircos[1]);

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
                TI   = ndof*ii;
                TIp1 = TI+1;

                fact1 = NN[ii]*Jmod;

                resi[TI]   += fact1*res(0);
                resi[TIp1] += fact1*res(1);

                for(jj=0;jj<nlbf;jj++)
                {
                  TJ   = ndof*jj;
                  TJp1 = TJ+1;

                  fact = fact1 * NN[jj];

                  stiffness_local[TI][TJ]      +=  fact;
                  stiffness_local[TIp1][TJp1]  +=  fact;
                }
              }
           }
           //cout << " side4 done " << endl;
        }
//printForceVector();
    }
//printStiffnessMatrix();
//printf("\n\n");
//printForceVector();
*/
  return 0;
}





int NurbsElem2DNavierStokes3dof::calcStiffnessMatrix(double dt)
{

  return 0;
}


int NurbsElem2DNavierStokes3dof::calcMassMatrix(int lumpInd, double dt)
{

  return 0;
}





int NurbsElem2DNavierStokes3dof::calcOutput(double u1, double v1)
{

  return 0;
}










