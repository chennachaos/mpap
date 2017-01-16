
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




inline void  tempfunc3(double aa, double bb, MatrixXd& b, MatrixXd& ff)
{
    double  b1 = 2.0 - (2.0+aa*(bb-1.0))/bb;
    double  b2 = (aa - 2.0)/bb;

    ff(0,0) = b1*b(0,0)+(aa/bb)*(b(1,1)+b(2,2));
    ff(0,1) = (2.0-2.0/bb)*b(0,1);
    ff(0,2) = (1.0-aa)*b(0,1);
    ff(0,3) = b(1,1);

    ff(0,4) = (-2.0/bb)*b(0,1);
    ff(0,5) = (-aa*(bb-1.0)/bb)*b(0,0)+b2*b(1,1)+(aa/bb)*b(2,2);
    ff(0,6) = b(0,0);
    ff(0,7) = (1.0-aa)*b(0,1);

    ff(1,0) = (1.0-aa)*b(0,1);
    ff(1,1) = b(1,1);
    ff(1,2) = b2*b(0,0)+(-aa*(bb-1.0)/bb)*b(1,1)+(aa/bb)*b(2,2);
    ff(1,3) = (-2.0/bb)*b(0,1);

    ff(1,4) = b(0,0);
    ff(1,5) = (1.0-aa)*b(0,1);
    ff(1,6) = (2.0-2.0/bb)*b(0,1);
    ff(1,7) = b1*b(1,1)+(aa/bb)*(b(0,0)+b(2,2));
}



inline void  tempfunc4(int nlbf, double J, double* dircos, double* NN, double* dN_dx, double* dN_dy, MatrixXd& ff, MatrixXd& D)
{
    int  ii, TI, TIp1, TIp2;

    double BULK, lambda, mu;

    BULK = 4.009417e5;

    BULK = 200.471;
    mu   = 92.525;

    lambda = BULK - 2.0*mu/3.0;
    //BULK = lambda+mu;

    for(ii=0;ii<nlbf;ii++)
    {
        TI   = 3*ii;
        TIp1 = TI+1;
        TIp2 = TI+2;

        D(TI,0)    = (ff(0,0)*dN_dx[ii] + ff(0,1)*dN_dy[ii])*dircos[0];
        D(TI,0)   += (ff(0,2)*dN_dx[ii] + ff(0,3)*dN_dy[ii])*dircos[1];
        D(TIp1,0)  = (ff(0,4)*dN_dx[ii] + ff(0,5)*dN_dy[ii])*dircos[0];
        D(TIp1,0) += (ff(0,6)*dN_dx[ii] + ff(0,7)*dN_dy[ii])*dircos[1];
        D(TIp2,0)  = NN[ii]*dircos[0];

        D(TI,1)    = (ff(1,0)*dN_dx[ii] + ff(1,1)*dN_dy[ii])*dircos[0];
        D(TI,1)   += (ff(1,2)*dN_dx[ii] + ff(1,3)*dN_dy[ii])*dircos[1];
        D(TIp1,1)  = (ff(1,4)*dN_dx[ii] + ff(1,5)*dN_dy[ii])*dircos[0];
        D(TIp1,1) += (ff(1,6)*dN_dx[ii] + ff(1,7)*dN_dy[ii])*dircos[1];
        D(TIp2,1)  = NN[ii]*dircos[1];

        D(TI,2)   =  0.5*(J+1.0/J)*dN_dx[ii];
        D(TIp1,2) =  0.5*(J+1.0/J)*dN_dy[ii];
        D(TIp2,2) = -NN[ii]/BULK;

        //D(TI,2)   = 0.0;
        //D(TIp1,2) = 0.0;
        //D(TIp2,2) = 0.0;
    }

  return;
}


int NurbsElem2DStructSolidLSFEM3dof::calcStiffnessAndResidual2()
{
/*
   int ii, jj, gp1, gp2, TI, TIp1, TIp2, count1, index;

   double  fact, dvol, dvol0, Jac, dt, lambda, BULK, mu, pres, J, fact1, fact2, aa, bb, rho;
   double  b1, b2, b3, ci, N[nlbf], dN_dx[nlbf], dN_dy[nlbf];
   double  d2N_dx2[nlbf], d2N_dy2[nlbf], d2N_dxy[nlbf], d2N_dyx[nlbf];

   VectorXd  DJ(2), f(3), dp(3);
   MatrixXd  D(nsize,3), F(3,3), G(3,3), DF(5,2), b(3,3), Db(4,2), bdev(3,3), II(3,3);
   MatrixXd  D2u(2,3), ff(2,10), gg(2,4);

   II.setZero();
   II(0,0) = II(1,1) = II(2,2) = 1.0;

   bool pout = false;
   pout = (bool) matDat[3];

   BULK = matDat[0];
   mu = matDat[1];

   lambda = BULK - 2.0*mu/3.0;
   //BULK = lambda+mu;

   for(ii=0;ii<nsize;ii++)
     stiffness_local[ii].zero();
   
   resi.zero();

   double  *gaussweights = &(surf0->gaussweights[0]);
   double  *values1 = &(surf1->Values[0][0]);
   double  *values2 = &(surf1->Values[1][0]);
   double  *values3 = &(surf1->Values[2][0]);
   
   int *tt = &(surf0->IEN[elenum][0]);

   aa = matDat[4];
   bb = matDat[5];

   count1 = 0;
   for(gp2=0;gp2<nGP2;gp2++)
   {
   for(gp1=0;gp1<nGP1;gp1++)
   {
        index = count1*2;
        surf0->ShapeFunDerivatives2(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, d2N_dx2, d2N_dy2, d2N_dxy, d2N_dyx, Jac);

        D2u.setZero();
        F.setZero();
        for(ii=0;ii<nlbf;ii++)
        {
           TI = tt[ii];
             
           b1 = values1[TI];
           b2 = values2[TI];
           
           F(0,0) += b1*dN_dx[ii];
           F(0,1) += b1*dN_dy[ii];
           F(1,0) += b2*dN_dx[ii];
           F(1,1) += b2*dN_dy[ii];

           D2u(0,0) += b1*d2N_dx2[ii];
           D2u(0,1) += b1*d2N_dxy[ii];
           D2u(0,2) += b1*d2N_dy2[ii];

           D2u(1,0) += b2*d2N_dx2[ii];
           D2u(1,1) += b2*d2N_dxy[ii];
           D2u(1,2) += b2*d2N_dy2[ii];
        }

        F(0,0) += 1.0;
        F(1,1) += 1.0;
        F(2,2) = 1.0;

        b = F*F.transpose();
        G = F.inverse();

        if(CompareDoubles(bb, 3.0))
          b(2,2) = 1.0;
        if(CompareDoubles(bb, 2.0))
          b(2,2) = 0.0;

        DF(0,0) = D2u(0,0)*G(0,0) + D2u(0,1)*G(1,0);
        DF(0,1) = D2u(0,0)*G(0,1) + D2u(0,1)*G(1,1);

        DF(1,0) = D2u(0,1)*G(0,0) + D2u(0,2)*G(1,0);
        DF(1,1) = D2u(0,1)*G(0,1) + D2u(0,2)*G(1,1);

        DF(2,0) = D2u(1,0)*G(0,0) + D2u(1,1)*G(1,0);
        DF(2,1) = D2u(1,0)*G(0,1) + D2u(1,1)*G(1,1);

        DF(3,0) = D2u(1,1)*G(0,0) + D2u(1,2)*G(1,0);
        DF(3,1) = D2u(1,1)*G(0,1) + D2u(1,2)*G(1,1);

        // derivatives of b matrix entries

        Db(0,0) = 2.0* (F(0,0)*DF(0,0) + F(0,1)*DF(1,0));
        Db(0,1) = 2.0* (F(0,0)*DF(0,1) + F(0,1)*DF(1,1));

        Db(1,0) = F(0,0)*DF(2,0) + F(0,1)*DF(3,0) + F(1,0)*DF(0,0) + F(1,1)*DF(1,0) ;
        Db(1,1) = F(0,0)*DF(2,1) + F(0,1)*DF(3,1) + F(1,0)*DF(0,1) + F(1,1)*DF(1,1) ;

        Db(2,0) = 2.0*(F(1,0)*DF(2,0) + F(1,1)*DF(3,0));
        Db(2,1) = 2.0*(F(1,0)*DF(2,1) + F(1,1)*DF(3,1));

        Db(3,0) = 0.0;
        Db(3,1) = 0.0;

        // gradient of J

        J = F(0,0)*F(1,1) - F(0,1)*F(1,0);

        DJ(0) = F(0,0)*DF(3,0) + F(1,1)*DF(0,0) - F(0,1)*DF(2,0) - F(1,0)*DF(1,0);
        DJ(1) = F(0,0)*DF(3,1) + F(1,1)*DF(0,1) - F(0,1)*DF(2,1) - F(1,0)*DF(1,1);

        fact1 = mu/pow(J,aa);
        fact2 = -aa*fact1/J;

        b1 = 2.0 - (2.0+aa*(bb-1.0))/bb;
        b2 = (aa - 2.0)/bb;

        ff(0,0) = b1*Db(0,0)+(aa/bb)*(Db(2,0)+Db(3,0)) + (1.0-aa)*Db(1,1);
        ff(0,1) = (2.0-2.0/bb)*Db(1,0) + Db(2,1);
        ff(0,2) = b1*b(0,0)+(aa/bb)*(b(1,1)+b(2,2));
        ff(0,3) = (3.0-aa-2.0/bb)*b(0,1);
        ff(0,4) = b(1,1);

        ff(0,5) = (-2.0/bb)*Db(1,0) + Db(0,1);
        ff(0,6) = (1.0-aa)*Db(1,1) - (aa*(bb-1.0)/bb)*Db(0,0)+b2*Db(2,0)+(aa/bb)*Db(3,0);
        ff(0,7) = (-2.0/bb)*b(0,1);
        ff(0,8) = (1.0-aa*(bb-1.0)/bb)*b(0,0) + b2*b(1,1)+(aa/bb)*b(2,2);
        ff(0,9) = (1.0-aa)*b(0,1);

        ff(1,0) = (1.0-aa)*Db(1,0)+b2*Db(0,1) - (aa*(bb-1.0)/bb)*Db(2,1)+(aa/bb)*Db(3,1);
        ff(1,1) = Db(2,0) + (-2.0/bb)*Db(1,1);
        ff(1,2) = (1.0-aa)*b(0,1);
        ff(1,3) = b2*b(0,0)+(1.0-aa*(bb-1.0)/bb)*b(1,1) + (aa/bb)*b(2,2);
        ff(1,4) = (-2.0/bb)*b(0,1);

        ff(1,5) = Db(0,0) + (2.0-2.0/bb)*Db(1,1);
        ff(1,6) = (1.0-aa)*Db(1,0)+(aa/bb)*(Db(0,1)+Db(3,1)) + b1*Db(2,1);
        ff(1,7) = b(0,0);
        ff(1,8) = (3.0-aa-2.0/bb)*b(0,1);
        ff(1,9) = (aa/bb)*(b(0,0)+b(2,2)) + b1*b(1,1);

        gg(0,0) = (b1*b(0,0)+(aa/bb)*(b(1,1)+b(2,2)))*DJ(0) + (1.0-aa)*b(0,1)*DJ(1);
        gg(0,1) = (2.0-2.0/bb)*b(0,1)*DJ(0) + b(1,1)*DJ(1);
        gg(0,2) = (-2.0*b(0,1)/bb)*DJ(0) + b(0,0)*DJ(1);
        gg(0,3) = ((-aa*(bb-1.0)/bb)*b(0,0)+b2*b(1,1)+(aa/bb)*b(2,2))*DJ(0) + (1.0-aa)*b(0,1)*DJ(1);

        gg(1,0) = (1.0-aa)*b(0,1)*DJ(0) + (b2*b(0,0)-(aa*(bb-1.0)/bb)*b(1,1)+(aa/bb)*b(2,2))*DJ(1);
        gg(1,1) = b(1,1)*DJ(0) + (-2.0/bb)*b(0,1)*DJ(1);
        gg(1,2) = b(0,0)*DJ(0) + (2.0-2.0/bb)*b(0,1)*DJ(1);
        gg(1,3) = (1.0-aa)*b(0,1)*DJ(0) + ((aa/bb)*(b(0,0)+b(2,2))+b1*b(1,1))*DJ(1);

        ff *= fact1;
        gg *= fact2;

        index = count1*2;
        surf1->ShapeFunDerivatives2(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, d2N_dx2, d2N_dy2, d2N_dxy, d2N_dyx, Jac);

        dvol = Jac * gaussweights[count1] * JacMultFact;
        count1++;

        pres = 0.0;
        dp.setZero();
        for(ii=0;ii<nlbf;ii++)
        {
           b3 = values3[tt[ii]];

           pres += (b3 * N[ii]);
           
           dp(0)  += (b3 * dN_dx[ii]);
           dp(1)  += (b3 * dN_dy[ii]);

           TI   =  3*ii;
           TIp1 =  TI+1;
           TIp2 =  TI+2;

           D(TI,0)    =  ff(0,0)*dN_dx[ii]+ff(0,1)*dN_dy[ii]+ff(0,2)*d2N_dx2[ii]+ff(0,3)*d2N_dxy[ii]+ff(0,4)*d2N_dy2[ii];
           D(TI,0)   += (gg(0,0)*dN_dx[ii]+gg(0,1)*dN_dy[ii]);
           D(TIp1,0)  =  ff(0,5)*dN_dx[ii]+ff(0,6)*dN_dy[ii]+ff(0,7)*d2N_dx2[ii]+ff(0,8)*d2N_dxy[ii]+ff(0,9)*d2N_dy2[ii];
           D(TIp1,0) += (gg(0,2)*dN_dx[ii]+gg(0,3)*dN_dy[ii]);
           D(TIp2,0)  = dN_dx[ii];

           D(TI,1)    =  ff(1,0)*dN_dx[ii]+ff(1,1)*dN_dy[ii]+ff(1,2)*d2N_dx2[ii]+ff(1,3)*d2N_dxy[ii]+ff(1,4)*d2N_dy2[ii];
           D(TI,1)   += (gg(1,0)*dN_dx[ii]+gg(1,1)*dN_dy[ii]);
           D(TIp1,1)  =  ff(1,5)*dN_dx[ii]+ff(1,6)*dN_dy[ii]+ff(1,7)*d2N_dx2[ii]+ff(1,8)*d2N_dxy[ii]+ff(1,9)*d2N_dy2[ii];
           D(TIp1,1) += (gg(1,2)*dN_dx[ii]+gg(1,3)*dN_dy[ii]);
           D(TIp2,1)  = dN_dy[ii];

           D(TI,2)   =  0.5*(J+1.0/J)*dN_dx[ii];
           D(TIp1,2) =  0.5*(J+1.0/J)*dN_dy[ii];
           D(TIp2,2) = -N[ii]/BULK;
        }

        rho = rho0*timeFunction[0].prop/J;

        f(0) = -bforce[0] * rho ;
        f(1) = -bforce[1] * rho ;
        //f(0) = f(1) = 0.0;
        f(2) = 0.0;

        bdev = b - (b.trace()/bb)*II;

        //f(0) -= (fact1*(((bb-1.0)*Db(0,0)-Db(2,0)-Db(3,0))/bb + Db(1,1)) + fact2*(((bb-1.0)*b(0,0)-b(1,1)-b(2,2))*DJ(0)/bb + b(0,1)*DJ(1)) + dp(0));
        //f(1) -= (fact1*(((bb-1.0)*Db(2,1)-Db(0,1)-Db(3,1))/bb + Db(1,0)) + fact2*(((bb-1.0)*b(1,1)-b(0,0)-b(2,2))*DJ(1)/bb + b(1,0)*DJ(0)) + dp(1));

        f(0) -= (fact1*(((bb-1.0)*Db(0,0)-Db(2,0)-Db(3,0))/bb + Db(1,1)) + fact2*(bdev(0,0)*DJ(0) + bdev(0,1)*DJ(1)) + dp(0));
        f(1) -= (fact1*(((bb-1.0)*Db(2,1)-Db(0,1)-Db(3,1))/bb + Db(1,0)) + fact2*(bdev(1,0)*DJ(0) + bdev(1,1)*DJ(1)) + dp(1));
        f(2) -= (0.5*(J-1.0/J)- pres/BULK);

        //printf("\t f vector  = %14.8f \t%14.8f \n", dp(0), dp(1));
        //printf("\t f vector  = %14.8f \t%14.8f \t%14.8f \n", f(0), f(1), f(2));

        for(ii=0;ii<nsize;ii++)
        {
           resi[ii] += dvol*(D(ii,0)*f(0)+D(ii,1)*f(1)+D(ii,2)*f(2));

           for(jj=0;jj<nsize;jj++)
             stiffness_local[ii][jj]  +=  dvol*(D(ii,0)*D(jj,0)+D(ii,1)*D(jj,1) + D(ii,2)*D(jj,2)) ;
        }
    }//gp1
    }//gp2

//printStiffnessMatrix();
//printf("\n\n");
if(pout)
{ printForceVector();
  printf("\n\n");
}
*/
   return 0;
}



/*
int NurbsElem2DStructSolidLSFEM3dof::calcLoadVector2()
{
   if(tracflag)
   {
      for(int ii=0;ii<nsize;ii++)
        stiffness_local[ii].zero();

      resi.zero();

      VectorXd  res(3), temp(3), Normal1(3), Normal2(3);
      MatrixXd  D(nsize, 3), F(3,3), b(3,3), stre(3,3), ff(2,8);
      Normal1.setZero();
      Normal2.setZero();
      temp.setZero();
    
      double  *values1 = &(surf1->Values[0][0]);
      double  *values2 = &(surf1->Values[1][0]);
      double  *values3 = &(surf1->Values[2][0]);
   
      double *gausspoints1 = &(surf0->gausspoints1[0]);
      double *gausspoints2 = &(surf0->gausspoints2[0]);
      double *gaussweights1 = &(surf0->gaussweights1[0]);
      double *gaussweights2 = &(surf0->gaussweights2[0]);

      int *tt = &(surf0->IEN[elenum][0]);

      ListArray<EPOINT>  SKL;

      EPOINT  EP;

      int  ngbf1 = surf0->ngbf1;
      int  ngbf2 = surf0->ngbf2;

      bool pout = false;

      pout = (bool) matDat[3];

      int p = surf0->p, q = surf0->q, ii, jj, gp, TI, TIp1, TIp2, index, TJ, TJp1, TJp2;
      double  NN[nlbf], dN_dx[nlbf], dN_dy[nlbf], params[2], dircos[2];
      double  Jac, fact, fact1, fact2, ALPHA, ALPHA1, lambda;
      double  J, J1, Jmod, BULK, mu, b1, b2, b3, val1, val2, pres, aa, bb;
        
      BULK = matDat[0];
      mu = matDat[1];
      lambda = BULK - 2.0*mu/3.0;
      //BULK = lambda+mu;


      ALPHA = matDat[2];
      ALPHA1 = matDat[3];

      ALPHA = ALPHA1 = 1.0;
      //ALPHA = 1000.0;
 
      aa = matDat[4];
      bb = matDat[5];

      //aa = 1.0;
      //bb = 3.0;

      // side #1
      if(!CompareDoubles(tracdata[0][0],7777) || !CompareDoubles(tracdata[0][1],7777))
      {
          double  N[p+1];

          ListArray<CPOINT>  Pw1;
          Pw1.setDim(ngbf1);
          for(ii=0;ii<ngbf1;ii++)
            Pw1[ii] = surf1->Pw[ii][0];

          NurbsCURVE   curve_temp(Pw1, surf1->U, p);

          val1 = 0.5*uvalues[2];
          val2 = 0.5*(uvalues[0]+uvalues[1]);

          params[1] = 0.0;
          for(gp=0;gp<nGP1;gp++)   // loop over Gauss points
          {
              params[0] = val1*gausspoints1[gp] + val2;

              NurbsShapeFunctions2DAlg2(surf0, startindex[0], startindex[1], gausspoints1[gp], -1.0, N, J1, dircos);

              Normal1(0) = dircos[0];
              Normal1(1) = dircos[1];

              b1 = tracdata[0][0] *timeFunction[0].prop;
              b2 = tracdata[0][1] *timeFunction[0].prop;

              res(0) = -b1 * Normal1(0) - b2 * Normal1(1);
              res(1) = -b1 * Normal1(1) + b2 * Normal1(0);

              curve_temp.CurveDerPointRat(params[0], 1, SKL);

              J1 = SKL[1].Norm() * val1;

              Jmod = J1 * gaussweights1[gp];

              EP = SKL[1]/SKL[1].Norm();

              Normal2(0) =  EP.y;
              Normal2(1) = -EP.x;

              surf0->ShapeFunDerivatives3(&(startindex[0]), params, NN, dN_dx, dN_dy, Jac);

              if(pout)
              {
                printf(" tracX and tracY ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n", res(0), res(1), J1, Jmod, Jac);
                printf(" dircos ... %12.6f \t %12.6f \n", dircos[0], dircos[1]);
              }

              surf1->deformationGradient2(tt, NN, dN_dx, dN_dy, F, b, pres);

              J = F(0,0)*F(1,1) - F(0,1)*F(1,0);

              if(CompareDoubles(bb, 3.0))
                b(2,2) = 1.0;
              if(CompareDoubles(bb, 2.0))
                b(2,2) = 0.0;

              surf1->ShapeFunDerivatives3(&(startindex[0]), params, NN, dN_dx, dN_dy, Jac);

              fact1 = mu/pow(J,aa);

              b *= fact1;

              tempfunc3(aa, bb, b, ff);
              tempfunc4(nlbf, J, &(Normal2(0)), NN, dN_dx, dN_dy, ff, D);

              stre = b;

              fact2 = -b.trace()/bb + pres;

              stre(0,0) += fact2;
              stre(1,1) += fact2;

              res(2) = 0.0;
              res(0) -= (stre(0,0)*Normal2(0) + stre(0,1)*Normal2(1));
              res(1) -= (stre(1,0)*Normal2(0) + stre(1,1)*Normal2(1));
              res(2) -= (0.5*(J-1.0/J)- pres/BULK);

              //printf(" F ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f \n", F(0,0), F(0,1), F(1,0), F(1,1));
              //printf(" b ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f \n", b(0,0), b(0,1), b(1,0), b(1,1));
              //printf(" stre ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f \n", stre(0,0), stre(0,1), stre(1,0), stre(1,1));
              //printf(" res ... %12.6f \t %12.6f \n", res(0), res(1));
              //printf(" Jmod ... %12.6f \n\n", Jmod);

              //printf(" %14.8f \t%14.8f \t%14.8f \t%14.8f \t%14.8f \t%14.8f \t%14.8f \n\n", F(0,0), F(0,1), F(1,0), F(1,1), pres, res(0), res(1));
              //printf(" res ... %12.6f \t %12.6f \n", res(0), res(1));

              Jmod *= ALPHA1;
              for(ii=0;ii<nsize;ii++)
              {
                resi[ii] += Jmod*(D(ii,0)*res(0)+D(ii,1)*res(1)+D(ii,2)*res(2));

                for(jj=0;jj<nsize;jj++)
                  stiffness_local[ii][jj]  +=  Jmod*(D(ii,0)*D(jj,0)+D(ii,1)*D(jj,1)+D(ii,2)*D(jj,2)) ;
              }
          }
          if(pout) cout << " side1 done " << endl;
        }
        // side #3
        if(!CompareDoubles(tracdata[2][0],7777) || !CompareDoubles(tracdata[2][1],7777))
        {
           double  N[p+1];

           ListArray<CPOINT>  Pw1;
           Pw1.setDim(ngbf1);
           for(ii=0;ii<ngbf1;ii++)
             Pw1[ii] = surf1->Pw[ii][ngbf2-1];

           NurbsCURVE   curve_temp(Pw1, surf1->U, p);

           val1 = 0.5*uvalues[2];
           val2 = 0.5*(uvalues[0]+uvalues[1]);

           params[1] = 1.0;
           for(gp=0;gp<nGP1;gp++)   // loop over Gauss points
           {
              params[0] = val1*gausspoints1[gp] + val2;

              NurbsShapeFunctions2DAlg2(surf0, startindex[0], startindex[1], gausspoints1[gp], 1.0, N, J1, dircos);

              dircos[0] *= -1.0;
              dircos[1] *= -1.0;

              Normal1(0) = dircos[0];
              Normal1(1) = dircos[1];

              b1 = tracdata[2][0] *timeFunction[0].prop;
              b2 = tracdata[2][1] *timeFunction[0].prop;

              res(0) = -b1 * Normal1(0) - b2 * Normal1(1);
              res(1) = -b1 * Normal1(1) + b2 * Normal1(0);

              curve_temp.CurveDerPointRat(params[0], 1, SKL);

              J1 = SKL[1].Norm() * val1;

              Jmod = J1 * gaussweights1[gp];

              EP = SKL[1]/SKL[1].Norm();

              Normal2(0) =  EP.y;
              Normal2(1) = -EP.x;

              surf0->ShapeFunDerivatives3(&(startindex[0]), params, NN, dN_dx, dN_dy, Jac);

              if(pout)
              {
                printf(" tracX and tracY ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n", res(0), res(1), J1, Jmod, Jac);
                printf(" dircos ... %12.6f \t %12.6f \n", dircos[0], dircos[1]);
              }

              surf1->deformationGradient2(tt, NN, dN_dx, dN_dy, F, b, pres);

              J = F(0,0)*F(1,1) - F(0,1)*F(1,0);

              if(CompareDoubles(bb, 3.0))
                b(2,2) = 1.0;
              if(CompareDoubles(bb, 2.0))
                b(2,2) = 0.0;

              surf1->ShapeFunDerivatives3(&(startindex[0]), params, NN, dN_dx, dN_dy, Jac);

              fact1 = mu/pow(J,aa);

              b *= fact1;

              tempfunc3(aa, bb, b, ff);
              tempfunc4(nlbf, J, &(Normal2(0)), NN, dN_dx, dN_dy, ff, D);

              stre = b;

              fact2 = -b.trace()/bb + pres;

              stre(0,0) += fact2;
              stre(1,1) += fact2;

              res(2) = 0.0;
              res(0) -= (stre(0,0)*Normal2(0) + stre(0,1)*Normal2(1));
              res(1) -= (stre(1,0)*Normal2(0) + stre(1,1)*Normal2(1));
              res(2) -= (0.5*(J-1.0/J)- pres/BULK);

              //printf(" F ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f \n", F(0,0), F(0,1), F(1,0), F(1,1));
              //printf(" b ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f \n", b(0,0), b(0,1), b(1,0), b(1,1));
              //printf(" stre ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f \n", stre(0,0), stre(0,1), stre(1,0), stre(1,1));
              //printf(" res ... %12.6f \t %12.6f \n", res(0), res(1));
              //printf(" Jmod ... %12.6f \n\n", Jmod);

              //printf(" %14.8f \t%14.8f \t%14.8f \t%14.8f \t%14.8f \t%14.8f \t%14.8f \n\n", F(0,0), F(0,1), F(1,0), F(1,1), pres, res(0), res(1));
              //printf(" res ... %12.6f \t %12.6f \n", res(0), res(1));

              Jmod *= ALPHA1;

              for(ii=0;ii<nsize;ii++)
              {
                resi[ii] += Jmod*(D(ii,0)*res(0)+D(ii,1)*res(1)+D(ii,2)*res(2));

                for(jj=0;jj<nsize;jj++)
                  stiffness_local[ii][jj]  +=  Jmod*(D(ii,0)*D(jj,0)+D(ii,1)*D(jj,1)+D(ii,2)*D(jj,2)) ;
              }
           }
           if(pout) cout << " side3 done " << endl;
        }
        // side #2
        if(!CompareDoubles(tracdata[1][0],7777) || !CompareDoubles(tracdata[1][1],7777))
        {
           double   N[q+1];

           ListArray<CPOINT>  Pw1;
           Pw1.setDim(ngbf2);
           for(ii=0;ii<ngbf2;ii++)
             Pw1[ii] = surf1->Pw[ngbf1-1][ii];

           NurbsCURVE   curve_temp(Pw1, surf1->V, q);

           val1 = 0.5*vvalues[2];
           val2 = 0.5*(vvalues[0]+vvalues[1]);

           params[0] = 1.0;
           for(gp=0;gp<nGP2;gp++)   // loop over Gauss points
           {
              params[1] = val1*gausspoints2[gp] + val2;

              NurbsShapeFunctions2DAlg2(surf0, startindex[0], startindex[1], 1.0, gausspoints2[gp], N, J1, dircos);

              Normal1(0) = dircos[0];
              Normal1(1) = dircos[1];

              b1 = tracdata[1][0] *timeFunction[0].prop;
              b2 = tracdata[1][1] *timeFunction[0].prop;

              res(0) = -b1 * Normal1(0) - b2 * Normal1(1);
              res(1) = -b1 * Normal1(1) + b2 * Normal1(0);

              curve_temp.CurveDerPointRat(params[1], 1, SKL);

              J1 = SKL[1].Norm() * val1;

              //printf(" J1 ... %12.6f \t %12.6f \n", J1, SKL[1].Norm());

              Jmod = J1 * gaussweights2[gp];

              EP = SKL[1]/SKL[1].Norm();

              Normal2(0) =  EP.y;
              Normal2(1) = -EP.x;

              surf0->ShapeFunDerivatives3(&(startindex[0]), params, NN, dN_dx, dN_dy, Jac);

              if(pout)
              {
                printf(" tracX and tracY ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n", res(0), res(1), J1, Jmod, Jac);
                printf(" dircos ... %12.6f \t %12.6f \n", dircos[0], dircos[1]);
              }

              surf1->deformationGradient2(tt, NN, dN_dx, dN_dy, F, b, pres);

              J = F(0,0)*F(1,1) - F(0,1)*F(1,0);

              if(CompareDoubles(bb, 3.0))
                b(2,2) = 1.0;
              if(CompareDoubles(bb, 2.0))
                b(2,2) = 0.0;

              surf1->ShapeFunDerivatives3(&(startindex[0]), params, NN, dN_dx, dN_dy, Jac);

              temp = b*Normal2;
              fact1 = Normal2.dot(temp);
              //fact1 = 1.0;
              res *= sqrt(fact1)/J;

              //temp = F*Normal1;
              //fact1 = Normal2.dot(temp);
              //res *= fact1/J;

              fact1 = mu/pow(J,aa);

              b *= fact1;

              tempfunc3(aa, bb, b, ff);
              tempfunc4(nlbf, J, &(Normal2(0)), NN, dN_dx, dN_dy, ff, D);

              stre = b;

              fact2 = -b.trace()/bb + pres;

              stre(0,0) += fact2;
              stre(1,1) += fact2;

              res(2) = 0.0;
              res(0) -= (stre(0,0)*Normal2(0) + stre(0,1)*Normal2(1));
              res(1) -= (stre(1,0)*Normal2(0) + stre(1,1)*Normal2(1));
              res(2) -= (0.5*(J-1.0/J)- pres/BULK);

              //printf(" F ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f \n", F(0,0), F(0,1), F(1,0), F(1,1));
              //printf(" b ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f \n", b(0,0), b(0,1), b(1,0), b(1,1));
              //printf(" stre ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f \n", stre(0,0), stre(0,1), stre(1,0), stre(1,1));
              //printf(" res ... %12.6f \t %12.6f \n", res(0), res(1));
              //printf(" Jmod ... %12.6f \n\n", Jmod);

              //printf(" res ... %12.6f \t %12.6f \n", res(0), res(1));

              Jmod *= ALPHA1;

              for(ii=0;ii<nsize;ii++)
              {
                resi[ii] += Jmod*(D(ii,0)*res(0)+D(ii,1)*res(1)+D(ii,2)*res(2));

                for(jj=0;jj<nsize;jj++)
                  stiffness_local[ii][jj]  +=  Jmod*(D(ii,0)*D(jj,0)+D(ii,1)*D(jj,1)+D(ii,2)*D(jj,2)) ;
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

              surf0->ShapeFunDerivatives3(&(startindex[0]), params, NN, dN_dx, dN_dy, Jac);

              Jmod = J * gaussweights2[gp] ;

              //printf(" tracX and tracY ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n", res(0), res(1), J, Jmod, Jac);
              //printf(" dircos ... %12.6f \t %12.6f \n", dircos[0], dircos[1]);

              surf1->deformationGradient2(tt, NN, dN_dx, dN_dy, F, b, pres);

              J = F(0,0)*F(1,1) - F(0,1)*F(1,0);

              surf1->ShapeFunDerivatives3(&(startindex[0]), params, NN, dN_dx, dN_dy, Jac);


              res(0) = tracdata[3][0];
              res(1) = tracdata[3][1];

              for(ii=0;ii<nlbf;ii++)
              {
                 index = tt[ii];

                 res(0) -= values1[index]*NN[ii];
                 res(1) -= values2[index]*NN[ii];
              }

              res(2) = 0.0;
              //res(2) -= (0.5*(J-1.0/J)- pres/BULK);

              Jmod *= ALPHA;

              for(ii=0;ii<nlbf;ii++)
              {
                  TI   = 3*ii;
                  TIp1 = TI+1;
                  TIp2 = TI+2;

                  D(TI,0)    = NN[ii];
                  D(TIp1,0)  = 0.0;
                  D(TIp2,0)  = 0.0;

                  D(TI,1)    = 0.0;
                  D(TIp1,1)  = NN[ii];
                  D(TIp2,1)  = 0.0;

                  //D(TI,2)   =  0.5*(J+1.0/J)*dN_dx[ii];
                  //D(TIp1,2) =  0.5*(J+1.0/J)*dN_dy[ii];
                  //D(TIp2,2) = -NN[ii]/BULK;
                  D(TI,2)   = 0.0;
                  D(TIp1,2) = 0.0;
                  D(TIp2,2) = 0.0;
              }

              for(ii=0;ii<nsize;ii++)
              {
                resi[ii] += Jmod*(D(ii,0)*res(0)+D(ii,1)*res(1)+D(ii,2)*res(2));

                for(jj=0;jj<nsize;jj++)
                  stiffness_local[ii][jj]  +=  Jmod*(D(ii,0)*D(jj,0)+D(ii,1)*D(jj,1)+D(ii,2)*D(jj,2)) ;
              }
           }
           if(pout) cout << " side4 done " << endl;
        }
//printForceVector();
    }
//printStiffnessMatrix();
//printf("\n\n");
//printForceVector();

  return 0;
}
*/



//
int NurbsElem2DStructSolidLSFEM3dof::calcLoadVector2()
{
/*
   if(tracflag)
   {
      for(int ii=0;ii<nsize;ii++)
        stiffness_local[ii].zero();

      resi.zero();

      VectorXd  res(2), Normal2(3), temp(3);
      MatrixXd  D(nsize, 2), F(3,3), b(3,3), stre(3,3), ff(2,8);
      Normal2.setZero();
    
      double  *values1 = &(surf1->Values[0][0]);
      double  *values2 = &(surf1->Values[1][0]);
      double  *values3 = &(surf1->Values[2][0]);
   
      double *gausspoints1 = &(surf0->gausspoints1[0]);
      double *gausspoints2 = &(surf0->gausspoints2[0]);
      double *gaussweights1 = &(surf0->gaussweights1[0]);
      double *gaussweights2 = &(surf0->gaussweights2[0]);

      int *tt = &(surf0->IEN[elenum][0]);

      ListArray<CPOINT>  Pw1;
      ListArray<EPOINT>  SKL;

      EPOINT  EP, Normal;

      int  ngbf1 = surf0->ngbf1;
      int  ngbf2 = surf0->ngbf2;

      bool pout = false;

      pout = (bool) matDat[3];

      int p = surf0->p, q = surf0->q, ii, jj, gp, TI, TIp1, TIp2, index, TJ, TJp1, TJp2;
      double  NN[nlbf], dN_dx[nlbf], dN_dy[nlbf], params[2], dircos[2];
      double  Jac, fact, fact1, fact2, ALPHA, ALPHA1;
      double  J, J1, Jmod, BULK, mu, b1, b2, b3, val1, val2, pres, aa, bb;
        
      BULK = matDat[0];
      mu = matDat[1];

      ALPHA = matDat[2];
      ALPHA1 = matDat[3];

      ALPHA = ALPHA1 = 1.0;
 
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

              NurbsShapeFunctions2DAlg2(surf0, startindex[0], startindex[1], gausspoints1[gp], -1.0, N, J1, dircos);

              res(0) = tracdata[0][0];
              res(1) = tracdata[0][1];

              Jmod = J1 * gaussweights1[gp];

              surf0->ShapeFunDerivatives(&(startindex[0]), params, NN, dN_dx, dN_dy, Jac);

              for(ii=0;ii<nlbf;ii++)
              {
                 index = tt[ii];

                 res(0) -= values1[index]*NN[ii];
                 res(1) -= values2[index]*NN[ii];
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
          if(pout) cout << " side1 done " << endl;
        }
        // side #3
        if(!CompareDoubles(tracdata[2][0],7777) || !CompareDoubles(tracdata[2][1],7777))
        {
           double  N[p+1];

           val1 = 0.5*uvalues[2];
           val2 = 0.5*(uvalues[0]+uvalues[1]);

           params[1] = 1.0;
           for(gp=0;gp<nGP1;gp++)   // loop over Gauss points
           {
              params[0] = val1*gausspoints1[gp] + val2;

              NurbsShapeFunctions2DAlg2(surf0, startindex[0], startindex[1], gausspoints1[gp], 1.0, N, J1, dircos);

              res(0) = tracdata[2][0];
              res(1) = tracdata[2][1];

              Jmod = J1 * gaussweights1[gp];

              surf0->ShapeFunDerivatives(&(startindex[0]), params, NN, dN_dx, dN_dy, Jac);

              for(ii=0;ii<nlbf;ii++)
              {
                 index = tt[ii];

                 res(0) -= values1[index]*NN[ii];
                 res(1) -= values2[index]*NN[ii];
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

              NurbsShapeFunctions2DAlg2(surf0, startindex[0], startindex[1], 1.0, gausspoints2[gp], N, J1, dircos);

              res(0) = tracdata[1][0];
              res(1) = tracdata[1][1];

              Jmod = J1 * gaussweights1[gp];

              surf0->ShapeFunDerivatives(&(startindex[0]), params, NN, dN_dx, dN_dy, Jac);

              for(ii=0;ii<nlbf;ii++)
              {
                 index = tt[ii];

                 res(0) -= values1[index]*NN[ii];
                 res(1) -= values2[index]*NN[ii];
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

              //printf(" tracX and tracY ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n", res(0), res(1), J, Jmod, Jac);
              //printf(" dircos ... %12.6f \t %12.6f \n", dircos[0], dircos[1]);

              res(0) = tracdata[3][0];
              res(1) = tracdata[3][1];

              for(ii=0;ii<nlbf;ii++)
              {
                 index = tt[ii];

                 res(0) -= values1[index]*NN[ii];
                 res(1) -= values2[index]*NN[ii];
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
           if(pout) cout << " side4 done " << endl;
        }
//printForceVector();
    }
//printStiffnessMatrix();
//printf("\n\n");
//printForceVector();
*/
  return 0;
}
//





/*
int NurbsElem2DStructSolidLSFEM3dof::calcLoadVector2()
{
   if(tracflag)
   {
      for(int ii=0;ii<nsize;ii++)
        stiffness_local[ii].zero();

      resi.zero();

      VectorXd  res(2), temp(3), Normal1(3), Normal2(3);
      MatrixXd  D(nsize, 2), F(3,3), b(3,3), stre(3,3), ff(2,8);
      Normal2.setZero();
    
      double  *values1 = &(surf1->Values[0][0]);
      double  *values2 = &(surf1->Values[1][0]);
      double  *values3 = &(surf1->Values[2][0]);
   
      double *gausspoints1 = &(surf0->gausspoints1[0]);
      double *gausspoints2 = &(surf0->gausspoints2[0]);
      double *gaussweights1 = &(surf0->gaussweights1[0]);
      double *gaussweights2 = &(surf0->gaussweights2[0]);

      int *tt = &(surf0->IEN[elenum][0]);

      ListArray<CPOINT>  Pw1;
      ListArray<EPOINT>  SKL;

      EPOINT  EP, Normal;

      int  ngbf1 = surf0->ngbf1;
      int  ngbf2 = surf0->ngbf2;

      bool pout = false;

      pout = (bool) matDat[3];

      int p = surf0->p, q = surf0->q, ii, jj, gp, TI, TIp1, TIp2, index, TJ, TJp1, TJp2;
      double  NN[nlbf], dN_dx[nlbf], dN_dy[nlbf], params[2], dircos[2];
      double  Jac, fact, fact1, fact2, ALPHA, ALPHA1;
      double  J, J1, Jmod, BULK, mu, b1, b2, b3, val1, val2, pres, aa, bb;
        
      BULK = matDat[0];
      mu = matDat[1];

      ALPHA = matDat[2];
      ALPHA1 = matDat[3];

      ALPHA = ALPHA1 = 1.0;
 
      aa = matDat[4];
      bb = matDat[5];

      //aa = 1.0;
      //bb = 3.0;

      // side #1
      if(!CompareDoubles(tracdata[0][0],7777) || !CompareDoubles(tracdata[0][1],7777))
      {
          double  N[p+1];

          Pw1.setDim(ngbf1);
          for(ii=0;ii<ngbf1;ii++)
            Pw1[ii] = surf1->Pw[ii][0];

          NurbsCURVE   curve_temp(Pw1, surf1->U, p);

          val1 = 0.5*uvalues[2];
          val2 = 0.5*(uvalues[0]+uvalues[1]);

          params[1] = 0.0;
          for(gp=0;gp<nGP1;gp++)   // loop over Gauss points
          {
              params[0] = val1*gausspoints1[gp] + val2;

              NurbsShapeFunctions2DAlg2(surf0, startindex[0], startindex[1], gausspoints1[gp], -1.0, N, J1, dircos);

              b1 = tracdata[0][0] *timeFunction[0].prop;
              b2 = tracdata[0][1] *timeFunction[0].prop;

              res(0) = b1 * (-dircos[0]) + b2 * (-dircos[1]);
              res(1) = b1 * (-dircos[1]) + b2 * (dircos[0]);

              curve_temp.CurveDerPointRat(params[0], 1, SKL);

              J1 = SKL[1].Norm() * val1;

              Jmod = J1 * gaussweights1[gp];

              EP = SKL[1]/SKL[1].Norm();

              Normal.x =  EP.y;
              Normal.y = -EP.x;

              dircos[0] = Normal.x;
              dircos[1] = Normal.y;

              surf0->ShapeFunDerivatives3(&(startindex[0]), params, NN, dN_dx, dN_dy, Jac);

              if(pout)
              {
                Normal.print2screen();
                printf(" tracX and tracY ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n", res(0), res(1), J1, Jmod, Jac);
                printf(" dircos ... %12.6f \t %12.6f \n", dircos[0], dircos[1]);
              }

              surf1->deformationGradient2(tt, NN, dN_dx, dN_dy, F, b, pres);

              J = F(0,0)*F(1,1) - F(0,1)*F(1,0);

              if(CompareDoubles(bb, 3.0))
                b(2,2) = 1.0;
              if(CompareDoubles(bb, 2.0))
                b(2,2) = 0.0;

              surf1->ShapeFunDerivatives3(&(startindex[0]), params, NN, dN_dx, dN_dy, Jac);

              fact1 = mu/pow(J,aa);

              b *= fact1;

              tempfunc3(aa, bb, b, ff);

              tempfunc4(nlbf, dircos, NN, dN_dx, dN_dy, ff, D);

              stre = b;

              fact2 = -b.trace()/bb + pres;

              stre(0,0) += fact2;
              stre(1,1) += fact2;
              
              res(0) -= (stre(0,0)*dircos[0] + stre(0,1)*dircos[1]);
              res(1) -= (stre(1,0)*dircos[0] + stre(1,1)*dircos[1]);

              printf(" F ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f \n", F(0,0), F(0,1), F(1,0), F(1,1));
              printf(" b ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f \n", b(0,0), b(0,1), b(1,0), b(1,1));
              printf(" stre ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f \n", stre(0,0), stre(0,1), stre(1,0), stre(1,1));
              printf(" res ... %12.6f \t %12.6f \n", res(0), res(1));
              printf(" Jmod ... %12.6f \n\n", Jmod);


              //printf(" %14.8f \t%14.8f \t%14.8f \t%14.8f \t%14.8f \t%14.8f \t%14.8f \n\n", F(0,0), F(0,1), F(1,0), F(1,1), pres, res(0), res(1));
              //printf(" res ... %12.6f \t %12.6f \n", res(0), res(1));

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
           double  N[p+1];

           Pw1.setDim(ngbf1);
           for(ii=0;ii<ngbf1;ii++)
             Pw1[ii] = surf1->Pw[ii][ngbf2-1];

           NurbsCURVE   curve_temp(Pw1, surf0->U, p);

           val1 = 0.5*uvalues[2];
           val2 = 0.5*(uvalues[0]+uvalues[1]);

           params[1] = 1.0;

           for(gp=0;gp<nGP1;gp++)   // loop over Gauss points
           {
              params[0] = val1*gausspoints1[gp] + val2;

              NurbsShapeFunctions2DAlg2(surf0, startindex[0], startindex[1], gausspoints1[gp], 1.0, N, J1, dircos);

              dircos[0] *= -1.0;
              dircos[1] *= -1.0;

              b1 = tracdata[2][0] *timeFunction[0].prop;
              b2 = tracdata[2][1] *timeFunction[0].prop;

              res(0) = b1 * (-dircos[0]) + b2 * (-dircos[1]);
              res(1) = b1 * (-dircos[1]) + b2 * (dircos[0]);

              curve_temp.CurveDerPointRat(params[0], 1, SKL);

              J1 = SKL[1].Norm() * val1;

              Jmod = J1 * gaussweights1[gp];

              EP = SKL[1]/SKL[1].Norm();

              Normal.x =  EP.y;
              Normal.y = -EP.x;

              dircos[0] = -Normal.x;
              dircos[1] = -Normal.y;

              surf0->ShapeFunDerivatives3(&(startindex[0]), params, NN, dN_dx, dN_dy, Jac);

              if(pout)
              {
                Normal.print2screen();
                printf(" tracX and tracY ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n", res(0), res(1), J1, Jmod, Jac);
                printf(" dircos ... %12.6f \t %12.6f \n", dircos[0], dircos[1]);
              }

              surf1->deformationGradient2(tt, NN, dN_dx, dN_dy, F, b, pres);

              J = F(0,0)*F(1,1) - F(0,1)*F(1,0);

              if(CompareDoubles(bb, 3.0))
                b(2,2) = 1.0;
              if(CompareDoubles(bb, 2.0))
                b(2,2) = 0.0;

              surf1->ShapeFunDerivatives3(&(startindex[0]), params, NN, dN_dx, dN_dy, Jac);

              fact1 = mu/pow(J,aa);

              b *= fact1;

              tempfunc3(aa, bb, b, ff);

              tempfunc4(nlbf, dircos, NN, dN_dx, dN_dy, ff, D);

              stre = b;

              fact2 = -b.trace()/bb + pres;

              stre(0,0) += fact2;
              stre(1,1) += fact2;

              res(0) -= (stre(0,0)*dircos[0] + stre(0,1)*dircos[1]);
              res(1) -= (stre(1,0)*dircos[0] + stre(1,1)*dircos[1]);

              printf(" F ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f \n", F(0,0), F(0,1), F(1,0), F(1,1));
              printf(" b ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f \n", b(0,0), b(0,1), b(1,0), b(1,1));
              printf(" stre ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f \n", stre(0,0), stre(0,1), stre(1,0), stre(1,1));
              printf(" res ... %12.6f \t %12.6f \n", res(0), res(1));
              printf(" Jmod ... %12.6f \n\n", Jmod);


              //printf(" %14.8f \t%14.8f \t%14.8f \t%14.8f \t%14.8f \t%14.8f \t%14.8f \n\n", F(0,0), F(0,1), F(1,0), F(1,1), pres, res(0), res(1));
              //printf(" res ... %12.6f \t %12.6f \n", res(0), res(1));

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

           Pw1.setDim(ngbf2);
           for(ii=0;ii<ngbf2;ii++)
             Pw1[ii] = surf1->Pw[ngbf1-1][ii];

           NurbsCURVE   curve_temp(Pw1, surf0->V, q);

           val1 = 0.5*vvalues[2];
           val2 = 0.5*(vvalues[0]+vvalues[1]);

           params[0] = 1.0;
           for(gp=0;gp<nGP2;gp++)   // loop over Gauss points
           {
              params[1] = val1*gausspoints2[gp] + val2;

              NurbsShapeFunctions2DAlg2(surf0, startindex[0], startindex[1], 1.0, gausspoints2[gp], N, J1, dircos);

              b1 = tracdata[1][0] *timeFunction[0].prop;
              b2 = tracdata[1][1] *timeFunction[0].prop;

              res(0) = b1 * (-dircos[0]) + b2 * (-dircos[1]);
              res(1) = b1 * (-dircos[1]) + b2 * (dircos[0]);

              curve_temp.CurveDerPointRat(params[1], 1, SKL);

              printf(" res ... %12.6f \t %12.6f \n", J1, SKL[1].Norm());

              J1 = SKL[1].Norm() * val1;

              Jmod = J1 * gaussweights2[gp];

              EP = SKL[1]/SKL[1].Norm();

              Normal.x =  EP.y;
              Normal.y = -EP.x;

              dircos[0] = Normal.x;
              dircos[1] = Normal.y;

              surf0->ShapeFunDerivatives3(&(startindex[0]), params, NN, dN_dx, dN_dy, Jac);

              if(pout)
              {
                Normal.print2screen();
                printf(" tracX and tracY ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n", res(0), res(1), J1, Jmod, Jac);
                printf(" dircos ... %12.6f \t %12.6f \n", dircos[0], dircos[1]);
              }

              surf1->deformationGradient2(tt, NN, dN_dx, dN_dy, F, b, pres);

              J = F(0,0)*F(1,1) - F(0,1)*F(1,0);

              if(CompareDoubles(bb, 3.0))
                b(2,2) = 1.0;
              if(CompareDoubles(bb, 2.0))
                b(2,2) = 0.0;

              surf1->ShapeFunDerivatives3(&(startindex[0]), params, NN, dN_dx, dN_dy, Jac);

              //temp = F*Normal1;
              //fact1 = Normal2.dot(temp);
              //fact1 = 1.0;
              //res *= fact1/J;

              fact1 = mu/pow(J,aa);

              b *= fact1;

              tempfunc3(aa, bb, b, ff);

              tempfunc4(nlbf, dircos, NN, dN_dx, dN_dy, ff, D);

              stre = b;

              fact2 = -b.trace()/bb + pres;

              stre(0,0) += fact2;
              stre(1,1) += fact2;

              res(0) -= (stre(0,0)*dircos[0] + stre(0,1)*dircos[1]);
              res(1) -= (stre(1,0)*dircos[0] + stre(1,1)*dircos[1]);

              printf(" F ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f \n", F(0,0), F(0,1), F(1,0), F(1,1));
              printf(" b ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f \n", b(0,0), b(0,1), b(1,0), b(1,1));
              printf(" stre ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f \n", stre(0,0), stre(0,1), stre(1,0), stre(1,1));
              printf(" res ... %12.6f \t %12.6f \n", res(0), res(1));
              printf(" Jmod ... %12.6f \n\n", Jmod);

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

              surf0->ShapeFunDerivatives3(&(startindex[0]), params, NN, dN_dx, dN_dy, Jac);

              Jmod = J * gaussweights2[gp] ;

              //printf(" tracX and tracY ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n", res(0), res(1), J, Jmod, Jac);
              //printf(" dircos ... %12.6f \t %12.6f \n", dircos[0], dircos[1]);

              res(0) = tracdata[3][0];
              res(1) = tracdata[3][1];

              for(ii=0;ii<nlbf;ii++)
              {
                 index = tt[ii];

                 res(0) -= values1[index]*NN[ii];
                 res(1) -= values2[index]*NN[ii];
              }

              Jmod *= ALPHA;

              for(ii=0;ii<nlbf;ii++)
              {
                TI   = 3*ii;
                TIp1 = TI+1;
                //TIp2 = TI+2;

                fact1 = NN[ii]*Jmod;

                resi[TI]   += fact1*res(0);
                resi[TIp1] += fact1*res(1);

                for(jj=0;jj<nlbf;jj++)
                {
                  TJ   = 3*jj;
                  TJp1 = TJ+1;
                  //TJp2 = TJ+2;

                  fact = fact1 * NN[jj];

                  stiffness_local[TI][TJ]      +=  fact;
                  stiffness_local[TIp1][TJp1]  +=  fact;
                }
              }           
           }
           if(pout) cout << " side4 done " << endl;
        }
//printForceVector();
    }
//printStiffnessMatrix();
//printf("\n\n");
//printForceVector();

  return 0;
}
*/



int NurbsElem2DStructSolidLSFEM3dof::calcInternalForces()
{
/*
   int ii, jj, gp1, gp2, TI, TIp1, TIp2, count1, index;

   double  fact, dvol, dvol0, Jac, dt, lambda, BULK, mu, pres, J, fact1, fact2, aa, bb;
   double  b1, b2, b3, ci, N[nlbf], dN_dx[nlbf], dN_dy[nlbf];
   double  d2N_dx2[nlbf], d2N_dy2[nlbf], d2N_dxy[nlbf], d2N_dyx[nlbf];

   VectorXd  DJ(2), f(3), dp(3);
   MatrixXd  D(nsize,3), F(3,3), G(2,2), DF(5,2), b(3,3), Db(4,2), bdev(3,3), II(3,3);
   MatrixXd  D2u(2,3), ff(2,10), gg(2,4);

   II.setZero();
   II(0,0) = II(1,1) = II(2,2) = 1.0;

   bool pout = false;
   pout = (bool) matDat[3];

   BULK = matDat[0];
   mu = matDat[1];

   resi.zero();

   double  *gaussweights = &(surf0->gaussweights[0]);
   double  *values1 = &(surf1->Values[0][0]);
   double  *values2 = &(surf1->Values[1][0]);
   double  *values3 = &(surf1->Values[2][0]);
   
   int *tt = &(surf0->IEN[elenum][0]);

   aa = matDat[4];
   bb = matDat[5];

   count1 = 0;
   for(gp2=0;gp2<nGP2;gp2++)
   {
   for(gp1=0;gp1<nGP1;gp1++)
   {
        index = count1*2;

        surf0->ShapeFunDerivatives2(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, d2N_dx2, d2N_dy2, d2N_dxy, d2N_dyx, Jac);

        //dvol0 = Jac * gaussweights[count1] * JacMultFact;

        D2u.setZero();
        F.setZero();
        for(ii=0;ii<nlbf;ii++)
        {
           TI = tt[ii];
             
           b1 = values1[TI];
           b2 = values2[TI];
           b3 = values3[TI];
           
           F(0,0) += b1*dN_dx[ii];
           F(0,1) += b1*dN_dy[ii];
           F(1,0) += b2*dN_dx[ii];
           F(1,1) += b2*dN_dy[ii];

           D2u(0,0) += b1*d2N_dx2[ii];
           D2u(0,1) += b1*d2N_dxy[ii];
           D2u(0,2) += b1*d2N_dy2[ii];

           D2u(1,0) += b2*d2N_dx2[ii];
           D2u(1,1) += b2*d2N_dxy[ii];
           D2u(1,2) += b2*d2N_dy2[ii];
        }

        F(0,0) += 1.0;
        F(1,1) += 1.0;
        F(2,2) = 1.0;

        b = F*F.transpose();
        G = F.inverse();

        if(CompareDoubles(bb, 3.0))
          b(2,2) = 1.0;
        if(CompareDoubles(bb, 2.0))
          b(2,2) = 0.0;

        DF(0,0) = D2u(0,0)*G(0,0) + D2u(0,1)*G(1,0);
        DF(0,1) = D2u(0,0)*G(0,1) + D2u(0,1)*G(1,1);

        DF(1,0) = D2u(0,1)*G(0,0) + D2u(0,2)*G(1,0);
        DF(1,1) = D2u(0,1)*G(0,1) + D2u(0,2)*G(1,1);

        DF(2,0) = D2u(1,0)*G(0,0) + D2u(1,1)*G(1,0);
        DF(2,1) = D2u(1,0)*G(0,1) + D2u(1,1)*G(1,1);

        DF(3,0) = D2u(1,1)*G(0,0) + D2u(1,2)*G(1,0);
        DF(3,1) = D2u(1,1)*G(0,1) + D2u(1,2)*G(1,1);

        // derivatives of b matrix entries

        Db(0,0) = 2.0* (F(0,0)*DF(0,0) + F(0,1)*DF(1,0));
        Db(0,1) = 2.0* (F(0,0)*DF(0,1) + F(0,1)*DF(1,1));

        Db(1,0) = F(0,0)*DF(2,0) + F(0,1)*DF(3,0) + F(1,0)*DF(0,0) + F(1,1)*DF(1,0) ;
        Db(1,1) = F(0,0)*DF(2,1) + F(0,1)*DF(3,1) + F(1,0)*DF(0,1) + F(1,1)*DF(1,1) ;

        Db(2,0) = 2.0* (F(1,0)*DF(2,0) + F(1,1)*DF(3,0));
        Db(2,1) = 2.0* (F(1,0)*DF(2,1) + F(1,1)*DF(3,1));

        Db(3,0) = 0.0;
        Db(3,1) = 0.0;

        // gradient of J

        J = F(0,0)*F(1,1) - F(0,1)*F(1,0);

        DJ(0) = F(0,0)*DF(3,0) + F(1,1)*DF(0,0) - F(0,1)*DF(2,0) - F(1,0)*DF(1,0);
        DJ(1) = F(0,0)*DF(3,1) + F(1,1)*DF(0,1) - F(0,1)*DF(2,1) - F(1,0)*DF(1,1);

        fact1 = mu/pow(J,aa);
        fact2 = -aa*fact1/J;

        b1 = 2.0 - (2.0+aa*(bb-1.0))/bb;
        b2 = (aa - 2.0)/bb;

        ff(0,0) = b1*Db(0,0)+(aa/bb)*(Db(2,0)+Db(3,0)) + (1.0-aa)*Db(1,1);
        ff(0,1) = (2.0-2.0/bb)*Db(1,0) + Db(2,1);
        ff(0,2) = b1*b(0,0)+(aa/bb)*(b(1,1)+b(2,2));
        ff(0,3) = (3.0-aa-2.0/bb)*b(0,1);
        ff(0,4) = b(1,1);

        ff(0,5) = (-2.0/bb)*Db(1,0) + Db(0,1);
        ff(0,6) = (1.0-aa)*Db(1,1) - (aa*(bb-1.0)/bb)*Db(0,0)+b2*Db(2,0)+(aa/bb)*Db(3,0);
        ff(0,7) = (-2.0/bb)*b(0,1);
        ff(0,8) = (1.0-aa*(bb-1.0)/bb)*b(0,0) + b2*b(1,1)+(aa/bb)*b(2,2);
        ff(0,9) = (1.0-aa)*b(0,1);

        ff(1,0) = (1.0-aa)*Db(1,0)+b2*Db(0,1) - (aa*(bb-1.0)/bb)*Db(2,1)+(aa/bb)*Db(3,1);
        ff(1,1) = Db(2,0) + (-2.0/bb)*Db(1,1);
        ff(1,2) = (1.0-aa)*b(0,1);
        ff(1,3) = b2*b(0,0)+(1.0-aa*(bb-1.0)/bb)*b(1,1) + (aa/bb)*b(2,2);
        ff(1,4) = (-2.0/bb)*b(0,1);

        ff(1,5) = Db(0,0) + (2.0-2.0/bb)*Db(1,1);
        ff(1,6) = (1.0-aa)*Db(1,0)+(aa/bb)*(Db(0,1)+Db(3,1)) + b1*Db(2,1);
        ff(1,7) = b(0,0);
        ff(1,8) = (3.0-aa-2.0/bb)*b(0,1);
        ff(1,9) = (aa/bb)*(b(0,0)+b(2,2)) + b1*b(1,1);

        gg(0,0) = (b1*b(0,0)+(aa/bb)*(b(1,1)+b(2,2)))*DJ(0) + (1.0-aa)*b(0,1)*DJ(1);
        gg(0,1) = (2.0-2.0/bb)*b(0,1)*DJ(0) + b(1,1)*DJ(1);
        gg(0,2) = (-2.0*b(0,1)/bb)*DJ(0) + b(0,0)*DJ(1);
        gg(0,3) = ((-aa*(bb-1.0)/bb)*b(0,0)+b2*b(1,1)+(aa/bb)*b(2,2))*DJ(0) + (1.0-aa)*b(0,1)*DJ(1);

        gg(1,0) = (1.0-aa)*b(0,1)*DJ(0) + (b2*b(0,0)-(aa*(bb-1.0)/bb)*b(1,1)+(aa/bb)*b(2,2))*DJ(1);
        gg(1,1) = b(1,1)*DJ(0) + (-2.0/bb)*b(0,1)*DJ(1);
        gg(1,2) = b(0,0)*DJ(0) + (2.0-2.0/bb)*b(0,1)*DJ(1);
        gg(1,3) = (1.0-aa)*b(0,1)*DJ(0) + ((aa/bb)*(b(0,0)+b(2,2))+b1*b(1,1))*DJ(1);

        ff *= fact1;
        gg *= fact2;

        index = count1*2;

        surf1->ShapeFunDerivatives2(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, d2N_dx2, d2N_dy2, d2N_dxy, d2N_dyx, Jac);

        dvol = Jac * gaussweights[count1] * JacMultFact;
        count1++;

        pres = 0.0;
        dp.setZero();
        for(ii=0;ii<nlbf;ii++)
        {
           b3 = values3[tt[ii]];

           pres += (b3 * N[ii]);
           
           dp(0)  += (b3 * dN_dx[ii]);
           dp(1)  += (b3 * dN_dy[ii]);

           TI   =  3*ii;
           TIp1 =  TI+1;
           TIp2 =  TI+2;

           D(TI,0)    =  ff(0,0)*dN_dx[ii]+ff(0,1)*dN_dy[ii]+ff(0,2)*d2N_dx2[ii]+ff(0,3)*d2N_dxy[ii]+ff(0,4)*d2N_dy2[ii];
           D(TI,0)   += (gg(0,0)*dN_dx[ii] + gg(0,1)*dN_dy[ii]);
           D(TIp1,0)  =  ff(0,5)*dN_dx[ii]+ff(0,6)*dN_dy[ii]+ff(0,7)*d2N_dx2[ii]+ff(0,8)*d2N_dxy[ii]+ff(0,9)*d2N_dy2[ii];
           D(TIp1,0) += (gg(0,2)*dN_dx[ii] + gg(0,3)*dN_dy[ii]);
           D(TIp2,0)  = dN_dx[ii];

           D(TI,1)    =  ff(1,0)*dN_dx[ii]+ff(1,1)*dN_dy[ii]+ff(1,2)*d2N_dx2[ii]+ff(1,3)*d2N_dxy[ii]+ff(1,4)*d2N_dy2[ii];
           D(TI,1)   += (gg(1,0)*dN_dx[ii] + gg(1,1)*dN_dy[ii]);
           D(TIp1,1)  =  ff(1,5)*dN_dx[ii]+ff(1,6)*dN_dy[ii]+ff(1,7)*d2N_dx2[ii]+ff(1,8)*d2N_dxy[ii]+ff(1,9)*d2N_dy2[ii];
           D(TIp1,1) += (gg(1,2)*dN_dx[ii] + gg(1,3)*dN_dy[ii]);
           D(TIp2,1)  = dN_dy[ii];

           D(TI,2)   =  0.5*(J+1.0/J)*dN_dx[ii];
           D(TIp1,2) =  0.5*(J+1.0/J)*dN_dy[ii];

           //D(TI,2)   =  J*dN_dx[ii];
           //D(TIp1,2) =  J*dN_dy[ii];

           D(TIp2,2) = -N[ii]/BULK;
        }

        //f(0) = -bforce[0] * rho0 ;
        //f(1) = -bforce[1] * rho0 ;
        f(0) = f(1) = 0.0;
        f(2) = 0.0;

        bdev = b - (b.trace()/bb)*II;

        //f(0) -= (fact1*(((bb-1.0)*Db(0,0)-Db(2,0)-Db(3,0))/bb + Db(1,1)) + fact2*(((bb-1.0)*b(0,0)-b(1,1)-b(2,2))*DJ(0)/bb + b(0,1)*DJ(1)) + dp(0));
        //f(1) -= (fact1*(((bb-1.0)*Db(2,1)-Db(0,1)-Db(3,1))/bb + Db(1,0)) + fact2*(((bb-1.0)*b(1,1)-b(0,0)-b(2,2))*DJ(1)/bb + b(1,0)*DJ(0)) + dp(1));

        f(0) -= (fact1*(((bb-1.0)*Db(0,0)-Db(2,0)-Db(3,0))/bb + Db(1,1)) + fact2*(bdev(0,0)*DJ(0) + bdev(0,1)*DJ(1)) + dp(0));
        f(1) -= (fact1*(((bb-1.0)*Db(2,1)-Db(0,1)-Db(3,1))/bb + Db(1,0)) + fact2*(bdev(1,0)*DJ(0) + bdev(1,1)*DJ(1)) + dp(1));
        f(2) -= (0.5*(J-1.0/J)- pres/BULK);
        //f(2) -= (J - 1.0 - pres/BULK);

        for(ii=0;ii<nsize;ii++)
           resi[ii] += dvol*(D(ii,0)*f(0)+D(ii,1)*f(1)+D(ii,2)*f(2));
    }//gp1
    }//gp2

//printStiffnessMatrix();
//printf("\n\n");
if(pout)
{ printForceVector();
  printf("\n\n");
}
   if(tracflag)
   {
      resi.zero();

      VectorXd  res(2);
      MatrixXd  D(nsize, 2), F(3,3), b(3,3), stre(3,3), ff(2,8);
    
      double  *values1 = &(surf1->Values[0][0]);
      double  *values2 = &(surf1->Values[1][0]);
      double  *values3 = &(surf1->Values[2][0]);
   
      double *gausspoints1 = &(surf0->gausspoints1[0]);
      double *gausspoints2 = &(surf0->gausspoints2[0]);
      double *gaussweights1 = &(surf0->gaussweights1[0]);
      double *gaussweights2 = &(surf0->gaussweights2[0]);

      int *tt = &(surf0->IEN[elenum][0]);

      ListArray<CPOINT>  Pw1;
      ListArray<EPOINT>  SKL;

      EPOINT  EP, Normal;

      int  ngbf1 = surf0->ngbf1;
      int  ngbf2 = surf0->ngbf2;

      bool pout = false;

      pout = (bool) matDat[3];

      int p = surf0->p, q = surf0->q, ii, jj, gp, TI, TIp1, TIp2, index, TJ, TJp1, TJp2;
      double  NN[nlbf], dN_dx[nlbf], dN_dy[nlbf], params[2], dircos[2];
      double  Jac, fact, fact1, fact2, ALPHA, ALPHA1;
      double  J, J1, Jmod, BULK, mu, b1, b2, b3, val1, val2, pres, aa, bb;
        
      BULK = matDat[0];
      mu = matDat[1];

      ALPHA = matDat[2];
      ALPHA1 = matDat[3];

      ALPHA = ALPHA1 = 1.0;
 
      aa = matDat[4];
      bb = matDat[5];

      // side #1
      if(!CompareDoubles(tracdata[0][0],7777) || !CompareDoubles(tracdata[0][1],7777))
      {
          double  N[p+1];

          Pw1.setDim(ngbf1);
          for(ii=0;ii<ngbf1;ii++)
            Pw1[ii] = surf1->Pw[ii][0];

          NurbsCURVE   curve_temp(Pw1, surf1->U, p);

          val1 = 0.5*uvalues[2];
          val2 = 0.5*(uvalues[0]+uvalues[1]);

          params[1] = 0.0;
          for(gp=0;gp<nGP1;gp++)   // loop over Gauss points
          {
              params[0] = val1*gausspoints1[gp] + val2;

              NurbsShapeFunctions2DAlg2(surf0, startindex[0], startindex[1], gausspoints1[gp], -1.0, N, J1, dircos);

              b1 = tracdata[0][0] *timeFunction[0].prop;
              b2 = tracdata[0][1] *timeFunction[0].prop;

              res(0) = b1 * (-dircos[0]) + b2 * (-dircos[1]);
              res(1) = b1 * (-dircos[1]) + b2 * (dircos[0]);

              curve_temp.CurveDerPointRat(params[0], 1, SKL);

              J1 = SKL[1].Norm() * val1;

              Jmod = J1 * gaussweights1[gp];

              EP = SKL[1]/SKL[1].Norm();

              Normal.x =  EP.y;
              Normal.y = -EP.x;

              dircos[0] = Normal.x;
              dircos[1] = Normal.y;

              surf0->ShapeFunDerivatives(&(startindex[0]), params, NN, dN_dx, dN_dy, Jac);

              if(pout)
              {
                Normal.print2screen();
                printf(" tracX and tracY ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n", res(0), res(1), J1, Jmod, Jac);
                printf(" dircos ... %12.6f \t %12.6f \n", dircos[0], dircos[1]);
              }

              surf1->deformationGradient2(tt, NN, dN_dx, dN_dy, F, b, pres);

              J = F(0,0)*F(1,1) - F(0,1)*F(1,0);

              if(CompareDoubles(bb, 3.0))
                b(2,2) = 1.0;
              if(CompareDoubles(bb, 2.0))
                b(2,2) = 0.0;

              surf1->ShapeFunDerivatives(&(startindex[0]), params, NN, dN_dx, dN_dy, Jac);

              fact1 = mu/pow(J,aa);

              b *= fact1;

              tempfunc3(aa, bb, b, ff);

              tempfunc4(nlbf, J, dircos, NN, dN_dx, dN_dy, ff, D);

              stre = b;

              fact2 = -b.trace()/bb + pres;

              stre(0,0) += fact2;
              stre(1,1) += fact2;
              
              res(0) -= (stre(0,0)*dircos[0] + stre(0,1)*dircos[1]);
              res(1) -= (stre(1,0)*dircos[0] + stre(1,1)*dircos[1]);

              //printf(" %14.8f \t%14.8f \t%14.8f \t%14.8f \t%14.8f \t%14.8f \t%14.8f \n\n", F(0,0), F(0,1), F(1,0), F(1,1), pres, res(0), res(1));
              //printf(" res ... %12.6f \t %12.6f \n", res(0), res(1));

              Jmod *= ALPHA1;
              for(ii=0;ii<nsize;ii++)
                resi[ii] += Jmod*(D(ii,0)*res(0)+D(ii,1)*res(1));
          }
          if(pout) cout << " side1 done " << endl;
        }
        // side #3
        if(!CompareDoubles(tracdata[2][0],7777) || !CompareDoubles(tracdata[2][1],7777))
        {
           double  N[p+1];

           Pw1.setDim(ngbf1);
           for(ii=0;ii<ngbf1;ii++)
             Pw1[ii] = surf1->Pw[ii][ngbf2-1];

           NurbsCURVE   curve_temp(Pw1, surf0->U, p);

           val1 = 0.5*uvalues[2];
           val2 = 0.5*(uvalues[0]+uvalues[1]);

           params[1] = 1.0;

           for(gp=0;gp<nGP1;gp++)   // loop over Gauss points
           {
              params[0] = val1*gausspoints1[gp] + val2;

              NurbsShapeFunctions2DAlg2(surf0, startindex[0], startindex[1], gausspoints1[gp], 1.0, N, J1, dircos);

              dircos[0] *= -1.0;
              dircos[1] *= -1.0;

              b1 = tracdata[2][0] *timeFunction[0].prop;
              b2 = tracdata[2][1] *timeFunction[0].prop;

              res(0) = b1 * (-dircos[0]) + b2 * (-dircos[1]);
              res(1) = b1 * (-dircos[1]) + b2 * (dircos[0]);

              curve_temp.CurveDerPointRat(params[0], 1, SKL);

              J1 = SKL[1].Norm() * val1;

              Jmod = J1 * gaussweights1[gp];

              EP = SKL[1]/SKL[1].Norm();

              Normal.x =  EP.y;
              Normal.y = -EP.x;

              dircos[0] = -Normal.x;
              dircos[1] = -Normal.y;

              surf0->ShapeFunDerivatives(&(startindex[0]), params, NN, dN_dx, dN_dy, Jac);

              if(pout)
              {
                Normal.print2screen();
                printf(" tracX and tracY ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n", res(0), res(1), J1, Jmod, Jac);
                printf(" dircos ... %12.6f \t %12.6f \n", dircos[0], dircos[1]);
              }

              surf1->deformationGradient2(tt, NN, dN_dx, dN_dy, F, b, pres);

              J = F(0,0)*F(1,1) - F(0,1)*F(1,0);

              if(CompareDoubles(bb, 3.0))
                b(2,2) = 1.0;
              if(CompareDoubles(bb, 2.0))
                b(2,2) = 0.0;

              surf1->ShapeFunDerivatives(&(startindex[0]), params, NN, dN_dx, dN_dy, Jac);

              fact1 = mu/pow(J,aa);

              b *= fact1;

              tempfunc3(aa, bb, b, ff);

              tempfunc4(nlbf, J, dircos, NN, dN_dx, dN_dy, ff, D);

              stre = b;

              fact2 = -b.trace()/bb + pres;

              stre(0,0) += fact2;
              stre(1,1) += fact2;

              res(0) -= (stre(0,0)*dircos[0] + stre(0,1)*dircos[1]);
              res(1) -= (stre(1,0)*dircos[0] + stre(1,1)*dircos[1]);

              //printf(" %14.8f \t%14.8f \t%14.8f \t%14.8f \t%14.8f \t%14.8f \t%14.8f \n\n", F(0,0), F(0,1), F(1,0), F(1,1), pres, res(0), res(1));
              //printf(" res ... %12.6f \t %12.6f \n", res(0), res(1));

              Jmod *= ALPHA1;

              for(ii=0;ii<nsize;ii++)
                resi[ii] += Jmod*(D(ii,0)*res(0)+D(ii,1)*res(1));
           }
           if(pout) cout << " side3 done " << endl;
        }
        // side #2
        if(!CompareDoubles(tracdata[1][0],7777) || !CompareDoubles(tracdata[1][1],7777))
        {
           double   N[q+1];

           Pw1.setDim(ngbf2);
           for(ii=0;ii<ngbf2;ii++)
             Pw1[ii] = surf1->Pw[ngbf1-1][ii];

           NurbsCURVE   curve_temp(Pw1, surf0->V, q);

           val1 = 0.5*vvalues[2];
           val2 = 0.5*(vvalues[0]+vvalues[1]);

           params[0] = 1.0;
           for(gp=0;gp<nGP2;gp++)   // loop over Gauss points
           {
              params[1] = val1*gausspoints2[gp] + val2;

              NurbsShapeFunctions2DAlg2(surf0, startindex[0], startindex[1], 1.0, gausspoints2[gp], N, J1, dircos);

              b1 = tracdata[1][0] *timeFunction[0].prop;
              b2 = tracdata[1][1] *timeFunction[0].prop;

              res(0) = b1 * (-dircos[0]) + b2 * (-dircos[1]);
              res(1) = b1 * (-dircos[1]) + b2 * (dircos[0]);

              curve_temp.CurveDerPointRat(params[1], 1, SKL);

              J1 = SKL[1].Norm() * val1;

              Jmod = J1 * gaussweights1[gp];

              EP = SKL[1]/SKL[1].Norm();

              Normal.x =  EP.y;
              Normal.y = -EP.x;

              dircos[0] = Normal.x;
              dircos[1] = Normal.y;

              surf0->ShapeFunDerivatives(&(startindex[0]), params, NN, dN_dx, dN_dy, Jac);

              if(pout)
              {
                Normal.print2screen();
                printf(" tracX and tracY ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n", res(0), res(1), J1, Jmod, Jac);
                printf(" dircos ... %12.6f \t %12.6f \n", dircos[0], dircos[1]);
              }

              surf1->deformationGradient2(tt, NN, dN_dx, dN_dy, F, b, pres);

              J = F(0,0)*F(1,1) - F(0,1)*F(1,0);

              if(CompareDoubles(bb, 3.0))
                b(2,2) = 1.0;
              if(CompareDoubles(bb, 2.0))
                b(2,2) = 0.0;

              surf1->ShapeFunDerivatives(&(startindex[0]), params, NN, dN_dx, dN_dy, Jac);

              fact1 = mu/pow(J,aa);

              b *= fact1;

              tempfunc3(aa, bb, b, ff);

              tempfunc4(nlbf, J, dircos, NN, dN_dx, dN_dy, ff, D);

              stre = b;

              fact2 = -b.trace()/bb + pres;

              stre(0,0) += fact2;
              stre(1,1) += fact2;

              //printf(" b.trace() = %12.6f \t %12.6f \n", b.trace(), stre(2,2));

              res(0) -= (stre(0,0)*dircos[0] + stre(0,1)*dircos[1]);
              res(1) -= (stre(1,0)*dircos[0] + stre(1,1)*dircos[1]);

              //printf(" res ... %12.6f \t %12.6f \n", res(0), res(1));

              Jmod *= ALPHA1;
              for(ii=0;ii<nsize;ii++)
                resi[ii] += Jmod*(D(ii,0)*res(0)+D(ii,1)*res(1));
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

              //dircos[0] *= -1.0;
              //dircos[1] *= -1.0;

              //res(0) = tracdata[3][0] * (-dircos[0]) + tracdata[3][1] * (-dircos[1]);
              //res(1) = tracdata[3][0] * (-dircos[1]) + tracdata[3][1] * (dircos[0]);

              //printf(" tracX and tracY ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n", res(0), res(1), J, Jmod, Jac);
              //printf(" dircos ... %12.6f \t %12.6f \n", dircos[0], dircos[1]);

              res(0) = tracdata[3][0];
              res(1) = tracdata[3][1];

              //pres = 0.0;
              for(ii=0;ii<nlbf;ii++)
              {
                 index = tt[ii];

                 res(0) -= values1[index]*NN[ii];
                 res(1) -= values2[index]*NN[ii];
                 //pres   -= values3[index]*NN[ii];
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
              }           
           }
           if(pout) cout << " side4 done " << endl;
        }
//printForceVector();
    }
//printStiffnessMatrix();
//printf("\n\n");
//printForceVector();
*/
   return 0;
}



