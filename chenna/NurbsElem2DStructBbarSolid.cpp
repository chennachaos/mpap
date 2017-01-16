
#include <math.h>
#include "Debug.h"
#include "Plot.h"
#include "FunctionsElement.h"
#include "MpapTime.h"
#include "NurbsElem2DStructBbarSolid.h"
#include "NurbsShapeFunctions.h"
#include <assert.h>

#include "ComputerTime.h"
#include <Eigen/Dense>


using namespace Eigen;
using namespace std;



extern Plot     plot;
extern MpapTime mpapTime;
extern ComputerTime       computerTime;



NurbsElem2DStructBbarSolid::NurbsElem2DStructBbarSolid()
{
  if (debug) cout << " constructor NurbsElem2DStructBbarSolid\n\n";

  // cout << " constructor NurbsElem2DStructBbarSolid\n\n";
}



NurbsElem2DStructBbarSolid::~NurbsElem2DStructBbarSolid()
{
  if (debug) cout << " destructor NurbsElem2DStructBbarSolid\n\n";

  // cout << " destructor NurbsElem2DStructBbarSolid\n\n";

}

int NurbsElem2DStructBbarSolid::calcStiffnessMatrix(double dt)
{

  return 0;
}



void NurbsElem2DStructBbarSolid::contourplot(int index, int nCol, double umin, double umax)
{
  return;
}

int NurbsElem2DStructBbarSolid::calcMassMatrix(int lumpInd, double dt)
{


  return 0;
}


inline void Bernstein_Basis(int deg, double u, double* NN)
{
	int ii, ind;
	for(ii=0;ii<=deg;ii++)
	{
		ind = deg-ii;
		NN[ii] = Bin(deg,ii)*pow(u,ii)*pow((1-u),ind);
	}

	return;
}


/*
void NurbsElem2DStructBbarSolid::matl02(bool flag, double* eps, double* stre, double** cc)
{
   double K, G, twoG, lambda, treps, press, epstr;

   K = matDat[0];
   G = matDat[1];

   twoG = 2 * G;
   lambda = K - twoG/3.0;

   epstr = (eps[0] + eps[1] + eps[2]);

   press = lambda*epstr;

   stre[0] = twoG*eps[0] + press;  //normal component sig11
   stre[1] = twoG*eps[1] + press;  //normal component sig22
   stre[3] = twoG*eps[3];		// shear component sig12
   stre[2] = twoG*eps[2] + press;  //normal component sig33

   for(int ii=0;ii<3;ii++)
   {
      for(int jj=0;jj<3;jj++)
         cc[ii][jj] = lambda;

      cc[ii][3] = 0.0;
      cc[3][ii] = 0.0;
      cc[ii][ii] += twoG;
   }
   cc[3][3] = G;

   return;
}
*/



void NurbsElem2DStructBbarSolid::matl02(bool flag, double* eps, double* stre, double** cc)
{
  double K, G, twoG, lambda, treps, press, epstr;

  int ii, jj;

  K = matDat[0];
  G = matDat[1];

  twoG = 2 * G;
  lambda = K - twoG/3.0;

  epstr = (eps[0] + eps[1] + eps[2]);

  treps = epstr;

  if(flag) treps = eps[4];

  press = K * treps - twoG*epstr/3.0;

//   press = lambda*epstr;

  stre[0] = twoG*eps[0] + press;
  stre[1] = twoG*eps[3];
  stre[3] = twoG*eps[1] + press;
  stre[2] = twoG*eps[2] + press;


  cc[0][0] =  4.0/3.0;   cc[0][1] = -2.0/3.0;   cc[0][2] = -2.0/3.0;    cc[0][3] = 0.0;
  cc[1][0] = -2.0/3.0;   cc[1][1] =  4.0/3.0;   cc[1][2] = -2.0/3.0;    cc[1][3] = 0.0;
  cc[2][0] = -2.0/3.0;   cc[2][1] = -2.0/3.0;   cc[2][2] =  4.0/3.0;    cc[2][3] = 0.0;
  cc[3][0] = 0.0;        cc[3][1] = 0.0;        cc[3][2] = 0.0;         cc[3][3] = 1.0;

  for(ii=0;ii<4;ii++)
  {
    for(jj=0;jj<4;jj++)
      cc[ii][jj] *= G;
  }

  return;
}




int NurbsElem2DStructBbarSolid::calcStiffnessAndResidual()
{
//   NurbsElem2DStructBbarSolid::calcStiffnessAndResidual1();

   NurbsElem2DStructBbarSolid::calcStiffnessAndResidual2();

//   NurbsElem2DStructBbarSolid::calcStiffnessAndResidual3();

  return 0;
}



int NurbsElem2DStructBbarSolid::calcStiffnessAndResidual1()
{
/*
  if(finite || axsy)
   {
     cout << " Large strains OR AXISYMMETRIC cases are not supported by NurbsElem2DStructBbarSolid  " << endl;
     return -1;
   }
   
   double K, G;

   K = matDat[0];
   G = matDat[1];   

   int  err, isw, count, count1, ll, ii, jj, twoI, twoJ, twoIp1, twoJp1, index, gp1, gp2;

   int  p = surf0->p,  q  = surf0->q,  sizeM = p*q;

   double F[4], Flow[4], detF, F33, fact, dvol0, dt, J, bb1, bb2, dummy;

   double cc[4][4], stre[4], bc[2][3], Bbar[nlbf][3][2], g[2][nlbf];

   double N[nGP][nlbf], dN_dx[nGP][nlbf], dN_dy[nGP][nlbf], volume[nGP], Nbar[nGP][sizeM];


   MatrixXd MAB(sizeM,sizeM), invMAB(sizeM,sizeM), Cmat1(sizeM,nlbf), Cmat2(sizeM,nlbf), Wmat1(sizeM,nlbf), Wmat2(sizeM,nlbf);

   MAB.setZero();
   invMAB.setZero();
   Cmat1.setZero();
   Cmat2.setZero();
   Wmat1.setZero();
   Wmat2.setZero();


   KNOTVECTOR  Ulow, Vlow;

   getLowerOrderKnotVector(surf1->U, surf1->p, Ulow);
   getLowerOrderKnotVector(surf1->V, surf1->q, Vlow);


   double *gaussweights = &(surf0->gaussweights[0]);

   count1 = 0;

   // loop over Gauss points
   for(gp2=0;gp2<nGP2;gp2++)
   {
      for(gp1=0;gp1<nGP1;gp1++)
      {
           index = count1*2;

           NurbsShapeFunctions2DAlg11(surf0, startindex[0], startindex[1], knotsAtGPs[index], knotsAtGPs[index+1], N[count1], dN_dx[count1], dN_dy[count1], J);

           BasisFuns2D(&(Ulow[0]), Ulow.n, p-1, &(Vlow[0]), Vlow.n, q-1, knotsAtGPs[index], knotsAtGPs[index+1], Nbar[count1]);

           dvol0 = J * gaussweights[count1] * thick;

           volume[count1] = dvol0;

           for(ii=0;ii<sizeM;ii++)
           {
               fact = Nbar[count1][ii] * dvol0;

               // compute Mass matrix matrix

               for(jj=0;jj<sizeM;jj++)
                  MAB(ii,jj) += fact * Nbar[count1][jj];

               // compute C matrix

               for(jj=0;jj<nlbf;jj++)
               {
                  Cmat1(ii,jj) += ( fact * dN_dx[count1][jj] );
                  Cmat2(ii,jj) += ( fact * dN_dy[count1][jj] );
               }
           }

           count1++;
       }
    }


        // Compute the inverse of MASS Matrix of the projected space

        invMAB = MAB.inverse();

        // compute W matrix W = invMAB * C'

        Wmat1 = invMAB * Cmat1;

        Wmat2 = invMAB * Cmat2;

    for(ii=0;ii<nsize;ii++)
      stiffness_local[ii].zero();

    resi.zero();

    count  = 1;
    count1 = 0;
    ll     = 0;
    err    = 0;
    isw    = 3;


   count1 = 0;

   // loop over Gauss points
   for(gp2=0;gp2<nGP2;gp2++)
   {
      for(gp1=0;gp1<nGP1;gp1++)
      {
         dvol0 = volume[count1];

         NurbsShapeFunctions2DAlg55(surf1, startindex[0], startindex[1], 1, dN_dx[count1], dN_dy[count1], F, detF);

         for(ii=0;ii<nlbf;ii++)
         {
            g[0][ii] = 0.0;
            g[1][ii] = 0.0;
            for(jj=0;jj<sizeM;jj++)
            {
               g[0][ii] += Wmat1(jj,ii) * Nbar[count1][jj];
               g[1][ii] += Wmat2(jj,ii) * Nbar[count1][jj];
            }
         }

         NurbsShapeFunctions2DAlg55(surf1, startindex[0], startindex[1], 1, g[0], g[1], Flow, detF);

         dummy = 0.5 * (Flow[0] + Flow[3] - F[0] - F[3]);

//        cout << " Pressure " << '\t' << 1000*0.5*(Flow[0]+Flow[3]) << '\t' << 1000*0.5*(F[0]+F[3]) << '\t' << 1000*dummy << endl;

         F[0] += dummy ;
         F[3] += dummy ;
         F33 = 1.0;

         dt = mpapTime.dt;
         matlib2d_(matDat, F, &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, (Element*) this);
         if(err !=0)           return 1;


//
   cc[0][0] =  4.0/3.0;
   cc[0][1] = -2.0/3.0;
   cc[0][2] = -2.0/3.0;   
   cc[1][0] = -2.0/3.0;
   cc[1][1] =  4.0/3.0;
   cc[1][2] = -2.0/3.0;
   cc[2][0] = -2.0/3.0;
   cc[2][1] = -2.0/3.0;
   cc[2][2] =  4.0/3.0;   
   cc[3][3] =  1.0;
   cc[0][3] = cc[1][3] = cc[2][3] = 0.0;
   cc[3][0] = cc[3][1] = cc[3][2] = 0.0;

   for(ii=0;ii<4;ii++)
   {
      for(jj=0;jj<4;jj++)
         cc[ii][jj] *= G;
   }
//
        for(ii=0;ii<4;ii++)
        {
          stre[ii] *= dvol0;
          for(jj=0;jj<4;jj++)
            cc[ii][jj] *= dvol0;
        }


        //==============================================
        // CALCULATE Bbar Matrices
        //==============================================


        for(ii=0;ii<nlbf;ii++)
        {
           bb1 = 0.5 * (g[0][ii] - dN_dx[count1][ii]);
           bb2 = 0.5 * (g[1][ii] - dN_dy[count1][ii]);
           
           //bb1 = bb2 = 0.0;
           
           Bbar[ii][0][0] = dN_dx[count1][ii] + bb1;
           Bbar[ii][1][0] = bb1;
           Bbar[ii][2][0] = dN_dy[count1][ii];

           Bbar[ii][0][1] = bb2;
           Bbar[ii][1][1] = dN_dy[count1][ii] + bb2;
           Bbar[ii][2][1] = dN_dx[count1][ii];
       }


        //==============================================
        // CALCULATE TANGENT STIFFNESS
        //==============================================

        for(ii=0;ii<nlbf;ii++)
        {
           bc[0][0] = (Bbar[ii][0][0] * cc[0][0] + Bbar[ii][1][0] * cc[1][0] + Bbar[ii][2][0] * cc[3][0]);
           bc[0][1] = (Bbar[ii][0][0] * cc[0][1] + Bbar[ii][1][0] * cc[1][1] + Bbar[ii][2][0] * cc[3][1]);
           bc[0][2] = (Bbar[ii][0][0] * cc[0][3] + Bbar[ii][1][0] * cc[1][3] + Bbar[ii][2][0] * cc[3][3]);

           bc[1][0] = (Bbar[ii][0][1] * cc[0][0] + Bbar[ii][1][1] * cc[1][0] + Bbar[ii][2][1] * cc[3][0]);
           bc[1][1] = (Bbar[ii][0][1] * cc[0][1] + Bbar[ii][1][1] * cc[1][1] + Bbar[ii][2][1] * cc[3][1]);
           bc[1][2] = (Bbar[ii][0][1] * cc[0][3] + Bbar[ii][1][1] * cc[1][3] + Bbar[ii][2][1] * cc[3][3]);

           twoI   = 2*ii;
           twoIp1 = twoI+1;

           for(jj=ii;jj<nlbf;jj++)
           {
              twoJ   = 2*jj;
              twoJp1 = twoJ+1;

              stiffness_local[twoI][twoJ]     +=  (bc[0][0] * Bbar[jj][0][0] + bc[0][1] * Bbar[jj][1][0] + bc[0][2] * Bbar[jj][2][0]) ;
              stiffness_local[twoI][twoJp1]   +=  (bc[0][0] * Bbar[jj][0][1] + bc[0][1] * Bbar[jj][1][1] + bc[0][2] * Bbar[jj][2][1]) ;
              stiffness_local[twoIp1][twoJ]   +=  (bc[1][0] * Bbar[jj][0][0] + bc[1][1] * Bbar[jj][1][0] + bc[1][2] * Bbar[jj][2][0]) ;
              stiffness_local[twoIp1][twoJp1] +=  (bc[1][0] * Bbar[jj][0][1] + bc[1][1] * Bbar[jj][1][1] + bc[1][2] * Bbar[jj][2][1]) ;

           }

           // compute the residual
//           resi[twoI]   -= (dN_dx[count1][ii]*stre[0] + dN_dy[count1][ii]*stre[3]) ;
//           resi[twoIp1] -= (dN_dx[count1][ii]*stre[3] + dN_dy[count1][ii]*stre[1]) ;

           resi[twoI]   -= (Bbar[ii][0][0]*stre[0] + Bbar[ii][1][0]*stre[1] + Bbar[ii][2][0]*stre[3]) ;
           resi[twoIp1] -= (Bbar[ii][0][1]*stre[0] + Bbar[ii][1][1]*stre[1] + Bbar[ii][2][1]*stre[3]) ;

        }

        count++;
        count1++;
        ll += nivGP;
    }
  }
  
  for(ii=0;ii<nsize;ii++)
  {
     for(jj=ii+1;jj<nsize;jj++)
     {
        stiffness_local[jj][ii] = stiffness_local[ii][jj];
     }
  }
*/
  return 0;
}





int NurbsElem2DStructBbarSolid::calcStiffnessAndResidual2()
{
/*
   if(finite || axsy)
   {
     cout << " Large strains OR AXISYMMETRIC cases are not supported by NurbsElem2DStructBbarSolid  " << endl;
     return -1;
   }
   
   double K, G;

   K = matDat[0];
   G = matDat[1];
   

   int   err, isw, count, count1, ll, ii, jj, twoI, twoJ, twoIp1, twoJp1, index, gp1, gp2;

   int  p  = surf0->p,  q  = surf0->q, sizeM = p*q;

   double   F[4], Flow[4], detF, F33, fact, dvol0, dt, J, bb1, bb2, dummy, totvol=0.0;

   double   cc[4][4], stre[4], bc[2][4], Bbar[nlbf][4][2], g[2][nlbf];

   double   volume[nGP], N[nGP][nlbf], dN_dx[nGP][nlbf], dN_dy[nGP][nlbf], Nbar[nGP][sizeM];


   MatrixXd MAB(sizeM,sizeM), invMAB(sizeM,sizeM), Cmat1(sizeM,nlbf), Cmat2(sizeM,nlbf), Wmat1(sizeM,nlbf), Wmat2(sizeM,nlbf);

   MAB.setZero();
   invMAB.setZero();
   Cmat1.setZero();
   Cmat2.setZero();
   Wmat1.setZero();
   Wmat2.setZero();

   KNOTVECTOR  Ulow, Vlow;

   getLowerOrderKnotVector(surf1->U, surf1->p, Ulow);
   getLowerOrderKnotVector(surf1->V, surf1->q, Vlow);


   double *gaussweights = &(surf0->gaussweights[0]);

   count1 = 0;
   // loop over Gauss points
   for(gp2=0;gp2<nGP2;gp2++)
   {
      for(gp1=0;gp1<nGP1;gp1++)
      {
           index = count1*2;

           NurbsShapeFunctions2DAlg11(surf0, startindex[0], startindex[1], knotsAtGPs[index], knotsAtGPs[index+1], N[count1], dN_dx[count1], dN_dy[count1], J);

           BasisFuns2D(&(Ulow[0]), Ulow.n, p-1, &(Vlow[0]), Vlow.n, q-1, knotsAtGPs[index], knotsAtGPs[index+1], Nbar[count1]);

           dvol0 = J * gaussweights[count1] * thick;
           
           totvol += dvol0;

           volume[count1] = dvol0;

           for(ii=0;ii<sizeM;ii++)
           {
               fact = Nbar[count1][ii] * dvol0;

               // compute Mass matrix matrix

               for(jj=0;jj<sizeM;jj++)
               {
                    MAB(ii,jj) += fact * Nbar[count1][jj];
               }

               // compute C matrix

               for(jj=0;jj<nlbf;jj++)
               {
                  Cmat1(ii,jj) += ( fact * dN_dx[count1][jj] );
                  Cmat2(ii,jj) += ( fact * dN_dy[count1][jj] );
               }
           }

           count1++;
       }
    }
    
//    cout << '\t' << " Total Volume for element # " << elenum << " is = " << totvol << endl; cout << endl;


        // Compute the inverse of MASS Matrix of the projected space

        //if(elenum == 0 )

        invMAB = MAB.inverse();

        // compute W matrix W = invMAB * C'

        Wmat1 = invMAB * Cmat1;

        Wmat2 = invMAB * Cmat2;

    for(ii=0;ii<nsize;ii++)
      stiffness_local[ii].zero();

    resi.zero();

    count  = 1;
    count1 = 0;
    ll     = 0;
    err    = 0;
    isw    = 3;


   count1 = 0;

   // loop over Gauss points
   for(gp2=0;gp2<nGP2;gp2++)
   {
      for(gp1=0;gp1<nGP1;gp1++)
      {
         dvol0 = volume[count1];

         NurbsShapeFunctions2DAlg55(surf1, startindex[0], startindex[1], 1, dN_dx[count1], dN_dy[count1], F, detF);

         for(ii=0;ii<nlbf;ii++)
         {
            g[0][ii] = 0.0;
            g[1][ii] = 0.0;
            for(jj=0;jj<sizeM;jj++)
            {
               g[0][ii] += Wmat1(jj,ii) * Nbar[count1][jj];
               g[1][ii] += Wmat2(jj,ii) * Nbar[count1][jj];
            }
         }

         NurbsShapeFunctions2DAlg55(surf1, startindex[0], startindex[1], 1, g[0], g[1], Flow, detF);

         dummy = (Flow[0] + Flow[3] - F[0] - F[3])/3.0;

         F[0] += dummy ;
         F[3] += dummy ;
         F33 = 1.0 + dummy ;


        dt = mpapTime.dt;
        matlib2d_(matDat, F, &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, (Element*) this);
        if(err !=0)           return 1;

//        printf(" stresses ");        printf("\t%12.8f\t%12.8f\t%12.8f\t%16.12f\n\n", stre[0], stre[1], stre[2], (stre[0]+stre[1]+stre[2])/3.0);

//
        cout << " Ev and pressure " << dummy << '\t' << (stre[0]+stre[1]+stre[2])/3.0 << endl;
     printf("\n");


   cc[0][0] =  4.0/3.0;
   cc[0][1] = -2.0/3.0;
   cc[0][2] = -2.0/3.0;   
   cc[1][0] = -2.0/3.0;
   cc[1][1] =  4.0/3.0;
   cc[1][2] = -2.0/3.0;
   cc[2][0] = -2.0/3.0;
   cc[2][1] = -2.0/3.0;
   cc[2][2] =  4.0/3.0;   
   cc[3][3] =  1.0;
   cc[0][3] = cc[1][3] = cc[2][3] = 0.0;
   cc[3][0] = cc[3][1] = cc[3][2] = 0.0;

   for(ii=0;ii<4;ii++)
   {
      for(jj=0;jj<4;jj++)
         cc[ii][jj] *= G;
   }
// 

        for(ii=0;ii<4;ii++)
        {
          stre[ii] *= dvol0;
          for(jj=0;jj<4;jj++)
            cc[ii][jj] *= dvol0;
        }

        //==============================================
        // CALCULATE Bbar Matrices
        //==============================================


        for(ii=0;ii<nlbf;ii++)
        {
           bb1 = (g[0][ii] - dN_dx[count1][ii])/3.0;
           bb2 = (g[1][ii] - dN_dy[count1][ii])/3.0;
           
          // bb1 = bb2 = 0.0;
           
           Bbar[ii][0][0] = dN_dx[count1][ii] + bb1;
           Bbar[ii][1][0] = bb1;
           Bbar[ii][2][0] = bb1;
           Bbar[ii][3][0] = dN_dy[count1][ii];

           Bbar[ii][0][1] = bb2;
           Bbar[ii][1][1] = dN_dy[count1][ii] + bb2;
           Bbar[ii][2][1] = bb2;
           Bbar[ii][3][1] = dN_dx[count1][ii];
       }


        //==============================================
        // CALCULATE TANGENT STIFFNESS
        //==============================================

        for(ii=0;ii<nlbf;ii++)
        {
           bc[0][0] = (Bbar[ii][0][0] * cc[0][0] + Bbar[ii][1][0] * cc[1][0] + Bbar[ii][2][0] * cc[2][0] + Bbar[ii][3][0] * cc[3][0]);
           bc[0][1] = (Bbar[ii][0][0] * cc[0][1] + Bbar[ii][1][0] * cc[1][1] + Bbar[ii][2][0] * cc[2][1] + Bbar[ii][3][0] * cc[3][1]);
           bc[0][2] = (Bbar[ii][0][0] * cc[0][2] + Bbar[ii][1][0] * cc[1][2] + Bbar[ii][2][0] * cc[2][2] + Bbar[ii][3][0] * cc[3][2]);
           bc[0][3] = (Bbar[ii][0][0] * cc[0][3] + Bbar[ii][1][0] * cc[1][3] + Bbar[ii][2][0] * cc[2][3] + Bbar[ii][3][0] * cc[3][3]);

           bc[1][0] = (Bbar[ii][0][1] * cc[0][0] + Bbar[ii][1][1] * cc[1][0] + Bbar[ii][2][1] * cc[2][0] + Bbar[ii][3][1] * cc[3][0]);
           bc[1][1] = (Bbar[ii][0][1] * cc[0][1] + Bbar[ii][1][1] * cc[1][1] + Bbar[ii][2][1] * cc[2][1] + Bbar[ii][3][1] * cc[3][1]);
           bc[1][2] = (Bbar[ii][0][1] * cc[0][2] + Bbar[ii][1][1] * cc[1][2] + Bbar[ii][2][1] * cc[2][2] + Bbar[ii][3][1] * cc[3][2]);
           bc[1][3] = (Bbar[ii][0][1] * cc[0][3] + Bbar[ii][1][1] * cc[1][3] + Bbar[ii][2][1] * cc[2][3] + Bbar[ii][3][1] * cc[3][3]);

           twoI   = 2*ii;
           twoIp1 = twoI+1;

           for(jj=0; jj<nlbf; jj++)
           {
              twoJ   = 2*jj;
              twoJp1 = twoJ+1;

              stiffness_local[twoI][twoJ]     +=  (bc[0][0] * Bbar[jj][0][0] + bc[0][1] * Bbar[jj][1][0] + bc[0][2] * Bbar[jj][2][0] + bc[0][3] * Bbar[jj][3][0]) ;
              stiffness_local[twoI][twoJp1]   +=  (bc[0][0] * Bbar[jj][0][1] + bc[0][1] * Bbar[jj][1][1] + bc[0][2] * Bbar[jj][2][1] + bc[0][3] * Bbar[jj][3][1]) ;
              stiffness_local[twoIp1][twoJ]   +=  (bc[1][0] * Bbar[jj][0][0] + bc[1][1] * Bbar[jj][1][0] + bc[1][2] * Bbar[jj][2][0] + bc[1][3] * Bbar[jj][3][0]) ;
              stiffness_local[twoIp1][twoJp1] +=  (bc[1][0] * Bbar[jj][0][1] + bc[1][1] * Bbar[jj][1][1] + bc[1][2] * Bbar[jj][2][1] + bc[1][3] * Bbar[jj][3][1]) ;

           }

           // compute the residual
//           resi[twoI]   -= (dN_dx[count1][ii]*stre[0] + dN_dy[count1][ii]*stre[3]) ;
//           resi[twoIp1] -= (dN_dx[count1][ii]*stre[3] + dN_dy[count1][ii]*stre[1]) ;

           resi[twoI]   -= (Bbar[ii][0][0]*stre[0] + Bbar[ii][1][0]*stre[1] + Bbar[ii][3][0]*stre[3] + Bbar[ii][2][0]*stre[2]) ;
           resi[twoIp1] -= (Bbar[ii][0][1]*stre[0] + Bbar[ii][1][1]*stre[1] + Bbar[ii][3][1]*stre[3] + Bbar[ii][2][1]*stre[2]) ;

        }

        count++;
        count1++;
        ll += nivGP;
    }
  }
*/
  return 0;
}



int NurbsElem2DStructBbarSolid::calcInternalForces()
{

  return 0;
}




int NurbsElem2DStructBbarSolid::calcOutput(double u1, double v1)
{

  return 0;
}







void NurbsElem2DStructBbarSolid::discreteContourplot(int vartype, int varindex, int index, int nCol, double umin, double umax)
{
  if(index > nivGP)
  {
     cout << '\t' << " Error in NurbsElem2DStructBbarSolid::contourplot " << endl;
     return;
  }

   double outval[100];

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

  du = (uvalues[1] - uvalues[0])/nGP1;
  dv = (vvalues[1] - vvalues[0])/nGP2;

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

  EPOINT *EP1, *EP2, *EP3, *EP4;

  double x1[2] = {0.0, 0.0}, x2[2] = {0.0, 0.0}, x3[2] = {0.0, 0.0}, x4[2] = {0.0, 0.0}, u1 = 0.0;

  int ind1, ind2;

  count=0;

  for(jj=0;jj<nGP2;jj++)
  {
      ind1 = (nGP1+1)*jj;
      ind2 = (nGP1+1)*(jj+1);
      
      for(ii=0;ii<nGP1;ii++)
      {
          u1 = outval[count];

          EP1 = &(S1[ind1+ii]);
          EP2 = &(S1[ind2+ii]);
          EP3 = &(S1[ind2+ii+1]);
          EP4 = &(S1[ind1+ii+1]);

          x1[0] = EP1->x; x1[1] = EP1->y;
          x2[0] = EP2->x; x2[1] = EP2->y;
          x3[0] = EP3->x; x3[1] = EP3->y;
          x4[0] = EP4->x; x4[1] = EP4->y;

          // contour plot for 1st triangle
          plot.triangleContourPlot(x1, x2, x3, u1, u1, u1, umin, umax, nCol);

          // contour plot for 2nd triangle
          plot.triangleContourPlot(x1, x3, x4, u1, u1, u1, umin, umax, nCol);

          count++;
      }
  }

  return;
}








void NurbsElem2DStructBbarSolid::projectToKnots(bool extrapolateFlag, int vartype, int varindex, int index)
{
/*
   vals2project[0] = intVar2[indx];
   vals2project[1] = intVar2[(nGP1-1)*nivGP+indx];
   vals2project[2] = intVar2[nGP1*(nGP2-1)*nivGP+indx];
   vals2project[3] = intVar2[(nGP1*nGP2-1)*nivGP+indx];
*/

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

//   cout << '\t' << vals2project << endl; cout << endl;

  return;
}




void NurbsElem2DStructBbarSolid::projectStress(int varindex, double* outval)
{
/*
   int nivEL = nGP * nivGP;
   for(int ii=0;ii<nivEL;ii++)
     intVar2[ii] = intVar1[ii];

    int   err, isw, count, count1, ll, ii, jj, twoI, twoJ, twoIp1, twoJp1, index, gp1, gp2;

   int  p  = surf0->p,  q  = surf0->q, sizeM = p*q;

   double   F[4], Flow[4], detF, F33, fact, dvol0, dt, J, bb1, bb2, dummy, totvol=0.0;

   double   cc[4][4], stre[4], bc[2][4], Bbar[nlbf][4][2], g[2][nlbf];

   double   volume[nGP], N[nGP][nlbf], dN_dx[nGP][nlbf], dN_dy[nGP][nlbf], Nbar[nGP][sizeM];


   MatrixXd MAB(sizeM,sizeM), invMAB(sizeM,sizeM), Cmat1(sizeM,nlbf), Cmat2(sizeM,nlbf), Wmat1(sizeM,nlbf), Wmat2(sizeM,nlbf);

   MAB.setZero();
   invMAB.setZero();
   Cmat1.setZero();
   Cmat2.setZero();
   Wmat1.setZero();
   Wmat2.setZero();


   KNOTVECTOR  Ulow, Vlow;

   getLowerOrderKnotVector(surf1->U, surf1->p, Ulow);
   getLowerOrderKnotVector(surf1->V, surf1->q, Vlow);

   double *gaussweights = &(surf0->gaussweights[0]);

   count1 = 0;

   // loop over Gauss points
   for(gp2=0;gp2<nGP2;gp2++)
   {
      for(gp1=0;gp1<nGP1;gp1++)
      {
           index = count1*2;

           NurbsShapeFunctions2DAlg11(surf0, startindex[0], startindex[1], knotsAtGPs[index], knotsAtGPs[index+1], N[count1], dN_dx[count1], dN_dy[count1], J);

           BasisFuns2D(&(Ulow[0]), Ulow.n, p-1, &(Vlow[0]), Vlow.n, q-1, knotsAtGPs[index], knotsAtGPs[index+1], Nbar[count1]);

           dvol0 = J * gaussweights[count1] * thick;
           
           totvol += dvol0;

           volume[count1] = dvol0;

           for(ii=0;ii<sizeM;ii++)
           {
               fact = Nbar[count1][ii] * dvol0;

               // compute Mass matrix matrix

               for(jj=0;jj<sizeM;jj++)
               {
                    MAB(ii,ii) += fact * Nbar[count1][jj];
               }

               // compute C matrix

               for(jj=0;jj<nlbf;jj++)
               {
                  Cmat1(ii,jj) += ( fact * dN_dx[count1][jj] );
                  Cmat2(ii,jj) += ( fact * dN_dy[count1][jj] );
               }
           }


           count1++;
       }
    }
    

        // Compute the inverse of MASS Matrix of the projected space


        invMAB = MAB.inverse();

        // compute W matrix W = invMAB * C'

        Wmat1 = invMAB * Cmat1;

        Wmat2 = invMAB * Cmat2;

    count  = 1;
    count1 = 0;
    ll     = 0;
    err    = 0;
    isw    = 3;


   count1 = 0;

   // loop over Gauss points
   for(gp2=0;gp2<nGP2;gp2++)
   {
      for(gp1=0;gp1<nGP1;gp1++)
      {
         dvol0 = volume[count1];

         NurbsShapeFunctions2DAlg55(surf1, startindex[0], startindex[1], 1, dN_dx[count1], dN_dy[count1], F, detF);

         for(ii=0;ii<nlbf;ii++)
         {
            g[0][ii] = 0.0;
            g[1][ii] = 0.0;
            for(jj=0;jj<sizeM;jj++)
            {
               g[0][ii] += Wmat1(jj,ii) * Nbar[count1][jj];
               g[1][ii] += Wmat2(jj,ii) * Nbar[count1][jj];
            }
         }

         NurbsShapeFunctions2DAlg55(surf1, startindex[0], startindex[1], 1, g[0], g[1], Flow, detF);

         dummy = (Flow[0] + Flow[3] - F[0] - F[3])/3.0;

         F[0] += dummy ;
         F[3] += dummy ;
         F33 = 1.0 + dummy ;



        dt = mpapTime.dt;
        matlib2d_(matDat, F, &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, (Element*) this);


        if(varindex < 4)
           outval[count1] = stre[varindex];
        else if(varindex == 4)
           outval[count1] = sqrt((pow(stre[0]-stre[1],2.0) + pow(stre[1]-stre[2], 2.0) + pow(stre[2]-stre[0], 2.0) + 6* stre[3]*stre[3])/2);
        else if(varindex == 5)
           outval[count1] = (stre[0]+stre[1]+stre[2])/3;
        else
        {
           cout << '\t' << "    NurbsElem2DStructBbarSolid::projectStress .... : Error in 'varindex' " << endl;
           return;
        }

        count++;
        count1++;
        ll += nivGP;

    }
  }
*/
  return;
}




void NurbsElem2DStructBbarSolid::projectStrain(int vartype, int varindex, double* outval)
{
/*
   int nivEL = nGP * nivGP;
   for(int ii=0;ii<nivEL;ii++)
     intVar2[ii] = intVar1[ii];

   bool  intVarFlag = (nivGP > 0);

   int   err, isw, count, count1, ll, ii, jj, kk, mm, twoI, twoJ, twoIp1, twoJp1, index, cnt, gp1, gp2;

   int  p  = surf0->p,
        q  = surf0->q;

   double F[4], Flow[4], detF=0.0, F33=0.0, fact, fact1, fact2, fact3, fact4, stre[4], bc[2][4];

   double cc[4][4], g[2][nlbf];

   double dvol, dvol0=0.0, utemp, vtemp, dt, J, volume[nGP];

   double N[nGP][nlbf], dN_dx[nGP][nlbf], dN_dy[nGP][nlbf];

   double udiff, vdiff, voltot, bb1, bb2, dummy, epsdilbar;

    udiff = uvalues[1] - uvalues[0];
    vdiff = vvalues[1] - vvalues[0];

    int  sizeM = p*q;

    double  Nbar[nGP][sizeM];


    MatrixXd MAB(sizeM,sizeM), invMAB(sizeM,sizeM), Cmat1(sizeM,nlbf), Cmat2(sizeM,nlbf), Wmat1(sizeM,nlbf), Wmat2(sizeM,nlbf);

    MAB.setZero();
    invMAB.setZero();
    Cmat1.setZero();
    Cmat2.setZero();
    Wmat1.setZero();
    Wmat2.setZero();

   KNOTVECTOR  Ulow, Vlow;

   getLowerOrderKnotVector(surf1->U, surf1->p, Ulow);
   getLowerOrderKnotVector(surf1->V, surf1->q, Vlow);


    for(ii=0;ii<4;ii++)
       for(jj=0;jj<4;jj++)
         cc[ii][jj] = 0.0;


    double *gaussweights = &(surf0->gaussweights[0]);

   count1 = 0;

   // loop over Gauss points
   for(gp2=0;gp2<nGP2;gp2++)
   {
      for(gp1=0;gp1<nGP1;gp1++)
      {
           index = count1*2;

           NurbsShapeFunctions2DAlg11(surf0, startindex[0], startindex[1], knotsAtGPs[index], knotsAtGPs[index+1], N[count1], dN_dx[count1], dN_dy[count1], J);

           BasisFuns2D(&(Ulow[0]), Ulow.n, p-1, &(Vlow[0]), Vlow.n, q-1, knotsAtGPs[index], knotsAtGPs[index+1], Nbar[count1]);

           dvol0 = J * gaussweights[count1] * thick;

           volume[count1] = dvol0;

           // compute Mass matrix matrix

           for(ii=0;ii<sizeM;ii++)
           {
               fact = Nbar[count1][ii] * dvol0;
               for(jj=0;jj<sizeM;jj++)
               {
                  MAB(ii,jj) += fact * Nbar[count1][jj];
               }
           }


           // compute C matrix

           for(ii=0;ii<sizeM;ii++)
           {
              fact = dvol0*Nbar[count1][ii];

              for(jj=0;jj<nlbf;jj++)
              {
                 Cmat1(ii,jj) += dN_dx[count1][jj] * fact;
                 Cmat2(ii,jj) += dN_dy[count1][jj] * fact;
              }
           }

           count1++;
       }
    }


        // Compute the inverse of MASS Matrix of the projected space

        invMAB = MAB.inverse();

        // compute W matrix W = invMAB * C'

        Wmat1 = invMAB * Cmat1;

        Wmat2 = invMAB * Cmat2;



    count  = 1;
    count1 = 0;
    ll     = 0;
    err    = 0;
    isw    = 3;

    count1 = 0;

   // loop over Gauss points
   for(gp2=0;gp2<nGP2;gp2++)
   {
      for(gp1=0;gp1<nGP1;gp1++)
      {
         index = count1*2;

         dvol0 = volume[count1];

         NurbsShapeFunctions2DAlg55(surf1, startindex[0], startindex[1], 1, dN_dx[count1], dN_dy[count1], F, detF);

         for(ii=0;ii<nlbf;ii++)
         {
            g[0][ii] = 0.0;
            g[1][ii] = 0.0;
            for(jj=0;jj<sizeM;jj++)
            {
               g[0][ii] += Wmat1(jj,ii) * Nbar[count1][jj];
               g[1][ii] += Wmat2(jj,ii) * Nbar[count1][jj];
            }
         }


         NurbsShapeFunctions2DAlg55(surf1, startindex[0], startindex[1], 1, g[0], g[1], Flow, detF);

         dummy = (Flow[0] + Flow[3] - F[0] - F[3])/3.0;

         F[0] += dummy ;
         F[3] += dummy ;
         F33 = 1.0 + dummy ;


        dt = mpapTime.dt;
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

    }
  }
*/
  return;
}




void NurbsElem2DStructBbarSolid::projectIntVar(int index, double* outval)
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







