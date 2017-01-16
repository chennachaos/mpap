
#include <math.h>
#include "Debug.h"
#include "MpapTime.h"
#include "Plot.h"
#include "NurbsElem2DStokes.h"
#include "NurbsShapeFunctions.h"
#include <assert.h>
#include "ComputerTime.h"
#include "Functions.h"

using namespace std;
using namespace Eigen;

extern ComputerTime       computerTime;
extern MpapTime mpapTime;
extern Plot plot;


NurbsElem2DStokes::NurbsElem2DStokes()
{
  if (debug) cout << " constructor NurbsElem2DStokes\n\n";
}



NurbsElem2DStokes::~NurbsElem2DStokes()
{
  if (debug) cout << " destructor NurbsElem2DStokes\n\n";
}



void NurbsElem2DStokes::createTractionDataVariable()
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



int NurbsElem2DStokes::calcStiffnessAndResidual()
{
  if(ndof == 3)
    calcStiffnessAndResidual1();
  else
    calcStiffnessAndResidual2();

  return 0;
}



int NurbsElem2DStokes::calcInternalForces()
{
/*
  if(ndof == 3)
    calcInternalForces1();
  else
    calcInternalForces2();
*/
  return 0;
}



/*
int NurbsElem2DStructMixed2field::calcStiffnessAndResidual1()
{
  // Stokes flow
  
  int  err, isw, count, count1, ll, ii, jj, twoI, twoJ, twoIp1, twoJp1, index, gp1, gp2, mm;

   int  sizep = surf2->nlbf;

   double  FF[4], detF, F33, fact, dvol, dt, Jac, dummy, pres, bb1, bb2, volstr;
   double  N[nlbf], dN_dx[nlbf], dN_dy[nlbf], Nbar[sizep];
   double  b1, b2, b3, b4, b5, rho, mu, tau, eps, trgradu;

   double *gaussweights = &(surf0->gaussweights[0]);

   VectorXd  res(2), dp(2), Du(2), vel(3), vel2(3), vectmp(nlbf), force(3);
   MatrixXd  F(2,2);

   rho =  elmDat[4];
   mu  =  elmDat[5];
   //tau =  uvalues[2] * vvalues[2] * elmDat[5]/(12.0*mu);
   //eps =  elmDat[6];

   //int   count1, ii, jj, gp1, gp2, index, TI, TIp1, TIp2, TIp3, TJ, TJp1, TJp2, TJp3;
   int *tt = &(surf0->IEN[elenum][0]);

   Klocal.setZero();
   Flocal.setZero();

   resi2.zero();

   Kup.resize(nsize, sizep);
   Kup.setZero();

   double  *values1 = &(surf1->Values[0][0]);
   double  *values2 = &(surf1->Values[1][0]);

   count = 1;   ll = 0;   err = 0;   isw = 3;
   dt = mpapTime.dt;

   count1 = 0;
   for(gp2=0;gp2<nGP2;gp2++)
   {
      for(gp1=0;gp1<nGP1;gp1++)
      {
          index = count1*2;

          surf0->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), dN_dx, dN_dy, Jac);

          dvol = Jac * gaussweights[count1]  * JacMultFact;

          surf0->ShapeFunctions(knotsAtGPs[index], knotsAtGPs[index+1], N);

          surf1->deformationGradient(startindex[0], startindex[1], 1, dN_dx, dN_dy, FF, detF);

          pres = surf2->computeValueAndShanpeFns(1, knotsAtGPs[index], knotsAtGPs[index+1], Nbar);

//        printf(" stresses ");        printf("\t%12.8f\t%12.8f\t%12.8f\t%12.8f\n\n", stre[0], stre[1], stre[2], pres);

          vel(0) = surf1->computeValue(1, knotsAtGPs[index], knotsAtGPs[index+1]);
          vel(1) = surf1->computeValue(2, knotsAtGPs[index], knotsAtGPs[index+1]);
        
          //res(0) = stokes.computeForce(0, knotsAtGPs[index], knotsAtGPs[index+1]);
          //res(1) = stokes.computeForce(1, knotsAtGPs[index], knotsAtGPs[index+1]);

          //printf(" \t %14.8f \t %14.8f \t %14.8f \t %14.8f \t %14.8f \t %14.8f\n", res(0), res(1), JacMultFact, Jac, fact, dvol0);

          F(0,0) = FF[0] - 1.0;
          F(0,1) = FF[2] ;
          F(1,0) = FF[1] ;
          F(1,1) = FF[3] - 1.0;
          
          trgradu = F.trace();

          force.setZero();
          //force(0) = analy.computeXForce(uu, vv);
          //force(1) = analy.computeYForce(uu, vv);

          //==============================================
          // CALCULATE TANGENT STIFFNESS and RESIDUAL
          //==============================================

          for(ii=0;ii<nlbf;ii++)
          {
             twoI   = 2*ii;
             twoIp1 = twoI+1;

             b1 = dN_dx[ii]*dvol;
             b2 = dN_dy[ii]*dvol;
             b3 = N[ii]*dvol;

             b4 = mu*b1;
             b5 = mu*b2;

             Flocal[twoI]   += (b3*force(0) - b4*F(0,0) - b5*F(0,1) + b1*pres);
             Flocal[twoPI]  += (b3*force(1) - b4*F(1,0) - b5*F(1,1) + b2*pres);

             for(jj=0;jj<nlbf;jj++)
             {
                twoJ   = 2*jj;
                twoJp1 = twoJ+1;

                fact = b4*dN_dx[jj] + b5*dN_dy[jj];

                Klocal(twoI,   twoJ)    +=  fact ;
                Klocal(twoIp1, twoJp1)  +=  fact ;
             }

             for(jj=0;jj<sizep;jj++)
             {
                Kup(twoI,   jj) -= ( b1 * Nbar[jj] );
                Kup(twoIp1, jj) -= ( b2 * Nbar[jj] );
             }
          }

          fact = trgradu * dvol;

          for(ii=0;ii<sizep;ii++)
             resi2[ii] += fact * Nbar[ii];

          count++;
          count1++;
          ll += nivGP;
      }
   }
   //printStiffnessMatrix();
   //cout << Kup << endl;
   //cout << "\n\n" << endl;

  return 0;
}
*/


/*
int NurbsElem2DStokes::calcStiffnessAndResidual1()
{
   Stokes2DEx1  stokes;

   double  fact, dvol0, Jac, rho, mu, b1, b2, b3, ci, cj, trgradu;

   double  N[nlbf], dN_dx[nlbf], d2N_dx2[nlbf], dN_dy[nlbf], d2N_dy2[nlbf];
   VectorXd  res(3), dp(3), Du(3);

   MatrixXd  D(nsize,3);

   rho =  elmDat[3];
   mu  =  elmDat[4];

   int   count1, ii, jj, gp1, gp2, index, TI, TIp1, TIp2, TIp3, TJ, TJp1, TJp2, TJp3;

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

        fact = gaussweights[count1] * JacMultFact;
        count1++;
        dvol0 = Jac * fact;

        res.setZero();
        dp.setZero();
        Du.setZero();
        trgradu = 0.0;
        
        //res(0) = res(1) = 0.0;

        //res(0) = stokes.computeForce(0, knotsAtGPs[index], knotsAtGPs[index+1]);
        //res(1) = stokes.computeForce(1, knotsAtGPs[index], knotsAtGPs[index+1]);

        //printf(" \t %14.8f \t %14.8f \t %14.8f \t %14.8f \t %14.8f \t %14.8f\n", res(0), res(1), JacMultFact, Jac, fact, dvol0);

        for(ii=0;ii<nlbf;ii++)
        {
           TI = tt[ii];
             
           b1 = values1[TI];
           b2 = values2[TI];
           b3 = values3[TI];

           Du(0)  += b1*(d2N_dx2[ii]+d2N_dy2[ii]);
           Du(1)  += b2*(d2N_dx2[ii]+d2N_dy2[ii]);
           
           trgradu +=  (b1*dN_dx[ii]+b2*dN_dy[ii]);
           
           dp(0)  += b3 * dN_dx[ii];
           dp(1)  += b3 * dN_dy[ii];

           /////////////////////////////////////

           TI   =  3*ii;
           TIp1 =  TI+1;
           TIp2 =  TI+2;

           b1 = dN_dx[ii];
           b2 = dN_dy[ii];
           b3 = N[ii];

           ci = - mu * (d2N_dx2[ii] + d2N_dy2[ii]);

           D(TI,0)   = ci;   D(TI,1)   = 0.0;   D(TI,2)   = b1;
           D(TIp1,0) = 0.0;  D(TIp1,1) = ci;    D(TIp1,2) = b2;
           D(TIp2,0) = b1;   D(TIp2,1) = b2;    D(TIp2,2) = 0.0;
        }

        res(0) += (-mu*Du(0) + dp(0));
        res(1) += (-mu*Du(1) + dp(1));
        res(2) += trgradu;

        for(ii=0;ii<nsize;ii++)
        {
           resi[ii] -= dvol0*(D(ii,0)*res(0)+D(ii,1)*res(1)+D(ii,2)*res(2));
           for(jj=0;jj<nsize;jj++)
             stiffness_local[ii][jj]  +=  dvol0*(D(ii,0)*D(jj,0)+D(ii,1)*D(jj,1) + D(ii,2)*D(jj,2)) ;
        }
  }//gp1
  }//gp2

  //printForceVector();
  //cout << endl;  cout << endl;

  return 0;
}
*/


int NurbsElem2DStokes::calcStiffnessAndResidual1()
{
/*
  Stokes2DEx1  stokes;

   double  fact, dvol, Jac, rho, mu, b1, b2, b3, b4, b5, ci, cj, trgradu, pres, eps, tau;
   double  N[nlbf], dN_dx[nlbf], d2N_dx2[nlbf], dN_dy[nlbf], d2N_dy2[nlbf], d2N_dxy[nlbf], d2N_dyx[nlbf];

   VectorXd  res(2), dp(2), Du(2), vel(3), vel2(3), vectmp(nlbf), force(3), stres(3);
   MatrixXd  D(nsize,2), F(2,2), FN(2,2);
   D.setZero();

   rho =  elmDat[3];
   mu  =  elmDat[4];
   tau =  uvalues[2] * vvalues[2] * elmDat[5]/(12.0*mu);
   eps =  elmDat[6];

   int   count1, ii, jj, gp1, gp2, index, TI, TIp1, TIp2, TIp3, TJ, TJp1, TJp2, TJp3;

   Klocal.setZero();
   Flocal.setZero();

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

        surf0->ShapeFunDerivatives2(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, d2N_dx2, d2N_dy2, d2N_dxy, d2N_dyx, Jac);

        fact = gaussweights[count1] * JacMultFact;
        count1++;
        dvol = Jac * fact;

        vel.setZero();
        res.setZero();
        dp.setZero();
        Du.setZero();
        F.setZero();
        trgradu = pres = 0.0;
        
        //res(0) = stokes.computeForce(0, knotsAtGPs[index], knotsAtGPs[index+1]);
        //res(1) = stokes.computeForce(1, knotsAtGPs[index], knotsAtGPs[index+1]);

        //printf(" \t %14.8f \t %14.8f \t %14.8f \t %14.8f \t %14.8f \t %14.8f\n", res(0), res(1), JacMultFact, Jac, fact, dvol0);

        for(ii=0;ii<nlbf;ii++)
        {
           TI = tt[ii];
             
           b1 = values1[TI];
           b2 = values2[TI];
           b3 = values3[TI];

           vel(0) += b1 * N[ii];
           vel(1) += b2 * N[ii];

           F(0,0) += b1 * dN_dx[ii];
           F(0,1) += b1 * dN_dy[ii];
           F(1,0) += b2 * dN_dx[ii];
           F(1,1) += b2 * dN_dy[ii];

           vectmp[ii] = d2N_dx2[ii] + d2N_dy2[ii];
           
           Du(0)  += b1*vectmp[ii];
           Du(1)  += b2*vectmp[ii];

           pres   += b3 * N[ii];
           dp(0)  += b3 * dN_dx[ii];
           dp(1)  += b3 * dN_dy[ii];
        }

        trgradu = F.trace();

        force.setZero();
        //force(0) = analy.computeXForce(uu, vv);
        //force(1) = analy.computeYForce(uu, vv);

          for(ii=0;ii<nlbf;ii++)
          {
             TI   = 3*ii;
             TIp1 = TI+1;
             TIp2 = TI+2;

             // GLS stabilisation term

             b1 = dN_dx[ii];
             b2 = dN_dy[ii];
             b3 = N[ii];

             fact = - mu*vectmp[ii];

             D(TI,  0) = fact ;
             D(TIp1,0) = 0.0;
             D(TIp2,0) = b1;

             D(TI,  1) = 0.0;
             D(TIp1,1) = fact ;
             D(TIp2,1) = b2;

             ////////////////////////////////////

             b1 = dN_dx[ii]*dvol;
             b2 = dN_dy[ii]*dvol;
             b3 = N[ii]*dvol;

             b4 = mu*b1;
             b5 = mu*b2;

             for(jj=0;jj<nlbf;jj++)
             {
               TJ   = 3*jj;
               TJp1 = TJ+1;
               TJp2 = TJ+2;

               // Diffusion term

               fact = b4*dN_dx[jj] + b5*dN_dy[jj];

               Klocal(TI,   TJ)   += fact;
               Klocal(TIp1, TJp1) += fact;

               // pressure term
               Klocal(TI,   TJp2) -= (b1*N[jj]);
               Klocal(TIp1, TJp2) -= (b2*N[jj]);

               // continuity term
               Klocal(TIp2, TJ)   -= (b3*dN_dx[jj]);
               Klocal(TIp2, TJp1) -= (b3*dN_dy[jj]);
               Klocal(TIp2, TJp2) += 0.0;
               //Klocal(TIp2, TJp2) += (b3*eps*N[jj]);
               //Klocal(TIp2, TJp2) += (tau*b1*dN_dx[jj]+tau*b2*dN_dy[jj]);
             }

             //Flocal(TI)   += (b3*force(0) - b1*stres(0) - b2*stres(2));
             //Flocal(TIp1) += (b3*force(1) - b1*stres(2) - b2*stres(1));
             Flocal(TI)   += (b3*force(0) - b4*F(0,0) - b5*F(0,1) + b1*pres);
             Flocal(TIp1) += (b3*force(1) - b4*F(1,0) - b5*F(1,1) + b2*pres);
             Flocal(TIp2) += (b3*trgradu);
             //Flocal(TIp2) -= ( tau*b1*dp(0) + tau*b2*dp(1));
          }

          //res = force;
          res.setZero();

          res(0) -= (-mu*Du(0) + dp(0));
          res(1) -= (-mu*Du(1) + dp(1));

          //res(0) -=  dp(0);
          //res(1) -=  dp(1);

          dvol *= tau;

          Klocal += (dvol*(D*D.transpose() ));
          Flocal += (D*(dvol*res));
  }//gp1
  }//gp2

  //cout << Klocal << endl;
  //printForceVector();
  //cout << endl;  cout << endl;
*/
  return 0;
}



int NurbsElem2DStokes::calcLoadVector()
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
	//ALPHA = 1.0e6;

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

              //res(0) = tracdata[3][0] *timeFunction[0].prop;
              //res(1) = tracdata[3][1] *timeFunction[0].prop;

              res(0) = kova1.computeValue(0, EP.x, EP.y);
              res(1) = kova1.computeValue(1, EP.x, EP.y);

              //res(0) = sin(2.0*PI*params[0])*cos(2.0*PI*params[1]);
              //res(1) = cos(2.0*PI*params[0])*sin(2.0*PI*params[1]);

              //printf(" tracX and tracY ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n", res(0), res(1), J, Jmod, params[1]);
              //printf(" dircos ... %12.6f \t %12.6f \n", dircos[0], dircos[1]);

              for(ii=0;ii<nlbf;ii++)
              {
                 //printf(" %5d \t %12.6f \n",ii, NN[ii]);
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

              res.setZero();

              //res(0) = tracdata[1][0] *timeFunction[0].prop;
              //res(1) = tracdata[1][1] *timeFunction[0].prop;

              res(0) = kova1.computeValue(0, EP.x, EP.y);
              res(1) = kova1.computeValue(1, EP.x, EP.y);

              //res(0) = sin(2.0*PI*params[0])*cos(2.0*PI*params[1]);
              //res(1) = cos(2.0*PI*params[0])*sin(2.0*PI*params[1]);

              //printf(" tracX and tracY ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n", res(0), res(1), J, Jmod, params[1]);
              //printf(" dircos ... %12.6f \t %12.6f \n", dircos[0], dircos[1]);

              for(ii=0;ii<nlbf;ii++)
              {
                 //printf(" %5d \t %12.6f \n",ii, NN[ii]);
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

              res.setZero();

              //res(0) = tracdata[0][0] *timeFunction[0].prop;
              //res(1) = tracdata[0][1] *timeFunction[0].prop;

              res(0) = kova1.computeValue(0, EP.x, EP.y);
              res(1) = kova1.computeValue(1, EP.x, EP.y);


              //res(0) = sin(2.0*PI*params[0])*cos(2.0*PI*params[1]);
              //res(1) = cos(2.0*PI*params[0])*sin(2.0*PI*params[1]);

              //printf(" tracX and tracY ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n", res(0), res(1), J, Jmod, params[0]);
              //printf(" dircos ... %12.6f \t %12.6f \n", dircos[0], dircos[1]);

              for(ii=0;ii<nlbf;ii++)
              {
                 //printf(" %5d \t %12.6f \n",ii, NN[ii]);
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

              res.setZero();

              //res(0) = tracdata[2][0] *timeFunction[0].prop;
              //res(1) = tracdata[2][1] *timeFunction[0].prop;

              res(0) = kova1.computeValue(0, EP.x, EP.y);
              res(1) = kova1.computeValue(1, EP.x, EP.y);


              //res(0) = sin(2.0*PI*params[0])*cos(2.0*PI*params[1]);
              //res(1) = cos(2.0*PI*params[0])*sin(2.0*PI*params[1]);

              //printf(" tracX and tracY ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n", res(0), res(1), J, Jmod, params[0]);
              //printf(" dircos ... %12.6f \t %12.6f \n", dircos[0], dircos[1]);

              for(ii=0;ii<nlbf;ii++)
              {
                 //printf(" %5d \t %12.6f \n",ii, NN[ii]);
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
//printForceVector();
    }
//printStiffnessMatrix();
//printf("\n\n");
//printForceVector();
*/
  return 0;
}

int NurbsElem2DStokes::calcStiffnessAndResidual2()
{
/*
  Stokes2DEx1  stokes;

   double  fact, dvol0, Jac, nu, s, b1, b2, b3, b4, w, ci;

   VectorXd  N(nlbf), dN_dx(nlbf), dN_dy(nlbf), vel(2), f(4), dp(2), dw(2);
   
   MatrixXd  F(2,2), D(nsize,4);

   nu = 1.0/elmDat[4];
   s  =  elmDat[5];

   int   count1, ii, jj, gp1, gp2, index, TI, TIp1, TIp2, TIp3;

   for(ii=0;ii<nsize;ii++)
     stiffness_local[ii].zero();

   resi.zero();

   double *gaussweights = &(surf0->gaussweights[0]);
   
   double  *values1 = &(surf1->Values[0][0]);
   double  *values2 = &(surf1->Values[1][0]);
   double  *values3 = &(surf1->Values[2][0]);
   double  *values4 = &(surf1->Values[3][0]);
   
   int *tt = &(surf0->IEN[elenum][0]);
   
   count1 = 0;
   for(gp2=0;gp2<nGP2;gp2++)
   {
   for(gp1=0;gp1<nGP1;gp1++)
   {
        index = count1*2;

        surf0->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), &N(0), &dN_dx(0), &dN_dy(0), Jac);

        fact = gaussweights[count1] * JacMultFact;
        count1++;
        dvol0 = Jac * fact;
        
        //for(ii=0;ii<nlbf;ii++)
          //printf(" \t %14.8f\t %14.8f\t%14.8f\n", N[ii], dN_dx[ii], dN_dy[ii]);
          
        f.setZero();
        vel.setZero();
        F.setZero();
        dp.setZero();
        dw.setZero();
        w=0.0;

        if( CompareDoubles(s, 0.0) )
          f(0) = f(1) = 0.0;
        else
        {
          f(0) = stokes.computeForce(0, knotsAtGPs[index], knotsAtGPs[index+1]);
          f(1) = stokes.computeForce(1, knotsAtGPs[index], knotsAtGPs[index+1]);
        }

        for(ii=0;ii<nlbf;ii++)
        {
           TI = tt[ii];
             
           b1 = values1[TI];
           b2 = values2[TI];
           b3 = values3[TI];
           b4 = values4[TI];
           
           vel(0) += b1 * N[ii];
           vel(1) += b2 * N[ii];
           w      += b4 * N[ii];
           
           dp(0)  += b3 * dN_dx[ii];
           dp(1)  += b3 * dN_dy[ii];

           dw(0)  += b4 * dN_dx[ii];
           dw(1)  += b4 * dN_dy[ii];

           F(0,0) += b1 * dN_dx[ii];
           F(0,1) += b1 * dN_dy[ii];
           F(1,0) += b2 * dN_dx[ii];
           F(1,1) += b2 * dN_dy[ii];

           /////////////////////////////////////

           TI   =  4*ii;
           TIp1 =  TI+1;
           TIp2 =  TI+2;
           TIp3 =  TI+3;

           b1 = dN_dx[ii];
           b2 = dN_dy[ii];
           b3 = N[ii];

           D(TI,0)    = 0.0;
           D(TIp1,0)  = 0.0;
           D(TIp2,0)  = b1;
           D(TIp3,0)  = nu*b2;

           D(TI,1)    = 0.0;
           D(TIp1,1)  = 0.0;
           D(TIp2,1)  = b2;
           D(TIp3,1)  = -nu*b1;
           
           D(TI,2)    = b1;
           D(TIp1,2)  = b2;
           D(TIp2,2)  = 0.0;
           D(TIp3,2)  = 0.0;

           D(TI,3)    =  b2;
           D(TIp1,3)  = -b1;
           D(TIp2,3)  = 0.0;
           D(TIp3,3)  = b3;
        }

        f(0) -= (nu*dw(1) + dp(0));
        f(1) -= (- nu*dw(0) + dp(1));
        f(2) -= (F(0,0) + F(1,1));
        f(3) -= (w + F(0,1) - F(1,0));

        for(ii=0;ii<nsize;ii++)
        {
           resi[ii] += dvol0*(D(ii,0)*f(0)+D(ii,1)*f(1)+D(ii,2)*f(2)+D(ii,3)*f(3));

           for(jj=0;jj<nsize;jj++)
             stiffness_local[ii][jj]  +=  dvol0*(D(ii,0)*D(jj,0)+D(ii,1)*D(jj,1) + D(ii,2)*D(jj,2) + D(ii,3)*D(jj,3)) ;
        }
  }//gp1
  }//gp2

  //printForceVector();  
  //cout << endl;  cout << endl;
*/
  return 0;
}



int NurbsElem2DStokes::calcStiffnessMatrix(double dt)
{
  return 0;
}


int NurbsElem2DStokes::calcMassMatrix(int lumpInd, double dt)
{

  return 0;
}






