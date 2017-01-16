
#include <math.h>
#include "Debug.h"
//#include "Plot.h"
//#include "FunctionsElement.h"
#include "MpapTime.h"
#include "NurbsElem2DStructFbarSolid.h"
#include "NurbsShapeFunctions.h"
#include <assert.h>
#include <Eigen/Dense>
#include "QuadratureUtil.h"

using namespace Eigen;
using namespace std;


//extern Plot     plot;
extern MpapTime mpapTime;



NurbsElem2DStructFbarSolid::NurbsElem2DStructFbarSolid()
{
  if (debug) cout << " constructor NurbsElem2DStructFbarSolid\n\n";

  // cout << " constructor NurbsElem2DStructFbarSolid\n\n";

}



NurbsElem2DStructFbarSolid::~NurbsElem2DStructFbarSolid()
{
  if (debug) cout << " destructor NurbsElem2DStructFbarSolid\n\n";

  // cout << " destructor NurbsElem2DStructFbarSolid\n\n";

}

int NurbsElem2DStructFbarSolid::calcStiffnessMatrix(double dt)
{

  return 0;
}



void NurbsElem2DStructFbarSolid::contourplot(int index, int nCol, double umin, double umax)
{
  return;
}

int NurbsElem2DStructFbarSolid::calcMassMatrix(int lumpInd, double dt)
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


int NurbsElem2DStructFbarSolid::calcStiffnessAndResidual()
{
//   NurbsElem2DStructFbarSolid::calcStiffnessAndResidual1();

//   NurbsElem2DStructFbarSolid::calcStiffnessAndResidual2();

   NurbsElem2DStructFbarSolid::calcStiffnessAndResidual3();

  return 0;
}





int NurbsElem2DStructFbarSolid::calcStiffnessAndResidual1()
{
/*
   if(!finite || axsy)
   {
     cout << " Small strains OR AXISYMMETRIC cases are not supported by NurbsElem2DStructFbarSolid  " << endl;
     return -1;
   }

   int   err, isw, count, count1, ll, ii, jj, kk, mm, twoI, twoJ, twoIp1, twoJp1, index, cnt, gp1, gp2;

   int  p  = surf0->p,
        q  = surf0->q;

   double F[nGP][4], Flow[4], detF[nGP], detF0, F33=0.0, fact, fact1, fact2, fact3, fact4, stre[4], bc[2][5], alpha;

   double cc[4][4], Bbar[nlbf][4][2], g[2][nlbf], gmN[2][nlbf];

   double dvol, dvol0=0.0, utemp, vtemp, dt, J, volume[nGP], fact5, fact6, fact7, fact8, fact9;

   double N[nGP][nlbf], dN_dx[nGP][nlbf], dN_dy[nGP][nlbf];

   double udiff, vdiff, bb1, bb2, dummy, r1d3 = 1.0/3.0;

   udiff = uvalues[1] - uvalues[0];
   vdiff = vvalues[1] - vvalues[0];

    // memory allocation for extra variables

    int  sizeM = p*q;

    double  Nbar[nGP][sizeM];


    MatrixXd  MAB(sizeM,sizeM), invMAB(sizeM,sizeM), Cmat1(sizeM,nlbf), Cmat2(sizeM,nlbf), Wmat1(sizeM,nlbf), Wmat2(sizeM,nlbf);

    double  detFbar, detJbar1[sizeM], detJbar2[sizeM], g2[nsize][nsize], g3[nsize][nsize];

    double  Cmat3[sizeM][nsize][nsize], Wmat3[sizeM][nsize][nsize], Cmat4[sizeM][nsize][nsize], Wmat4[sizeM][nsize][nsize]  ;

    MAB.setZero();
    invMAB.setZero();
    Cmat1.setZero();
    Cmat2.setZero();
    Wmat1.setZero();
    Wmat2.setZero();

   KNOTVECTOR  Ulow, Vlow;

   getLowerOrderKnotVector(surf1->U, surf1->p, Ulow);
   getLowerOrderKnotVector(surf1->V, surf1->q, Vlow);


    for(kk=0;kk<sizeM;kk++)
    {
       detJbar1[kk] = 0.0;
       detJbar2[kk] = 0.0;

       for(ii=0;ii<nsize;ii++)
       {
          for(jj=0;jj<nsize;jj++)
          {
             Cmat3[kk][ii][jj] = 0.0;
             Wmat3[kk][ii][jj] = 0.0;
             Cmat4[kk][ii][jj] = 0.0;
             Wmat4[kk][ii][jj] = 0.0;
          }
       }
    }

    for(ii=0;ii<4;ii++)
    {
       for(jj=0;jj<4;jj++)
       {
         cc[ii][jj] = 0.0;
       }
    }

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

           NurbsShapeFunctions2DAlg55(surf1, startindex[0], startindex[1], 1, dN_dx[count1], dN_dy[count1], F[count1], detF[count1]);

           NurbsShapeFunctions2DAlg11(surf1, startindex[0], startindex[1], knotsAtGPs[index], knotsAtGPs[index+1], N[count1], dN_dx[count1], dN_dy[count1], J);

           dvol0 = J * gaussweights[count1] * thick;

           volume[count1] = dvol0;


           dummy = pow(detF[count1],r1d3);

           // compute Mass matrix and (detF)^(1/3)

           for(ii=0;ii<sizeM;ii++)
           {
               fact = Nbar[count1][ii] * dvol0;

               detJbar1[ii] += fact * dummy;

               for(jj=0;jj<sizeM;jj++)
               {
                  MAB(ii,jj) += fact * Nbar[count1][jj];
               }
           }

           // compute C matrix

           for(ii=0;ii<sizeM;ii++)
           {
              fact = dvol0 * Nbar[count1][ii] * dummy;

              for(jj=0;jj<nlbf;jj++)
              {
                 Cmat1(ii,jj) += dN_dx[count1][jj] * fact;
                 Cmat2(ii,jj) += dN_dy[count1][jj] * fact;
              }
           }


           for(kk=0;kk<sizeM;kk++)
           {
              fact = dvol0 * Nbar[count1][kk] * dummy;

              for(ii=0;ii<nlbf;ii++)
              {
                 twoI   = 2*ii;
                 twoIp1 = twoI+1;

                 fact1 = fact * dN_dx[count1][ii];
                 fact2 = fact * dN_dy[count1][ii];

                 for(jj=0;jj<nlbf;jj++)
                 {
                    twoJ   = 2*jj;
                    twoJp1 = twoJ+1;

                    Cmat3[kk][twoI][twoJ]     +=  ( fact1 * dN_dx[count1][jj] );
                    Cmat3[kk][twoI][twoJp1]   +=  ( fact1 * dN_dy[count1][jj] );
                    Cmat3[kk][twoIp1][twoJ]   +=  ( fact2 * dN_dx[count1][jj] );
                    Cmat3[kk][twoIp1][twoJp1] +=  ( fact2 * dN_dy[count1][jj] );

                    Cmat4[kk][twoI][twoJ]     +=  ( fact1 * dN_dx[count1][jj] );
                    Cmat4[kk][twoI][twoJp1]   +=  ( fact2 * dN_dx[count1][jj] );
                    Cmat4[kk][twoIp1][twoJ]   +=  ( fact1 * dN_dy[count1][jj] );
                    Cmat4[kk][twoIp1][twoJp1] +=  ( fact2 * dN_dy[count1][jj] );

                 }
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
         {
            for(jj=0;jj<nsize;jj++)
            {
               for(kk=0;kk<sizeM;kk++)
               {
                  for(ll=0;ll<sizeM;ll++)
                  {
                     Wmat3[kk][ii][jj] += invMAB(kk,ll) * Cmat3[ll][ii][jj];
                     Wmat4[kk][ii][jj] += invMAB(kk,ll) * Cmat4[ll][ii][jj];
                  }
               }
            }
         }

         for(kk=0;kk<sizeM;kk++)
         {
            detJbar2[kk] = 0.0;
            for(ll=0;ll<sizeM;ll++)
               detJbar2[kk] += invMAB(kk,ll) * detJbar1[ll];
         }





    for(ii=0;ii<nsize;ii++)
      stiffness_local[ii].zero();

    resi.zero();

    count  = 1;
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

         detFbar = 0.0;
         for(ii=0;ii<sizeM;ii++)
            detFbar += detJbar2[ii] * Nbar[count1][ii];

         alpha = detFbar/pow(detF[count1], r1d3);


         for(ii=0;ii<nlbf;ii++)
         {
            g[0][ii] = 0.0;
            g[1][ii] = 0.0;
            for(jj=0;jj<sizeM;jj++)
            {
               g[0][ii] += Wmat1(jj,ii) * Nbar[count1][jj];
               g[1][ii] += Wmat2(jj,ii) * Nbar[count1][jj];
            }
            g[0][ii] /= detFbar;
            g[1][ii] /= detFbar;
         }


         for(ii=0;ii<nsize;ii++)
         {
            for(jj=0;jj<nsize;jj++)
            {
               g2[ii][jj] = 0.0;
               g3[ii][jj] = 0.0;

               for(kk=0;kk<sizeM;kk++)
               {
                   g2[ii][jj] +=  ( Nbar[count1][kk] * Wmat3[kk][ii][jj] );
                   g3[ii][jj] +=  ( Nbar[count1][kk] * Wmat4[kk][ii][jj] );
               }

               g2[ii][jj] /= detFbar;
               g3[ii][jj] /= detFbar;
            }
         }


         // calculate Fbar

         for(ii=0;ii<4;ii++)
             F[count1][ii] = alpha * F[count1][ii];

         F33 = alpha ;



        dt = mpapTime.dt;
        matlib2d_(matDat, F[count1], &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, (Element*) this);
        if(err !=0)           return 1;

        dvol0 = dvol0 * F33;


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
           gmN[0][ii] = g[0][ii] - dN_dx[count1][ii];
           gmN[1][ii] = g[1][ii] - dN_dy[count1][ii];

           bb1 = gmN[0][ii]/3.0;
           bb2 = gmN[1][ii]/3.0;

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

        // 1.) Material Contribution
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

           for(jj=0;jj<nlbf;jj++)
           {
              twoJ   = 2*jj;
              twoJp1 = twoJ+1;

              stiffness_local[twoI][twoJ]     +=  (bc[0][0] * Bbar[jj][0][0] + bc[0][1] * Bbar[jj][1][0] + bc[0][2] * Bbar[jj][2][0] + bc[0][3] * Bbar[jj][3][0]) ;
              stiffness_local[twoI][twoJp1]   +=  (bc[0][0] * Bbar[jj][0][1] + bc[0][1] * Bbar[jj][1][1] + bc[0][2] * Bbar[jj][2][1] + bc[0][3] * Bbar[jj][3][1]) ;
              stiffness_local[twoIp1][twoJ]   +=  (bc[1][0] * Bbar[jj][0][0] + bc[1][1] * Bbar[jj][1][0] + bc[1][2] * Bbar[jj][2][0] + bc[1][3] * Bbar[jj][3][0]) ;
              stiffness_local[twoIp1][twoJp1] +=  (bc[1][0] * Bbar[jj][0][1] + bc[1][1] * Bbar[jj][1][1] + bc[1][2] * Bbar[jj][2][1] + bc[1][3] * Bbar[jj][3][1]) ;

           }

//           resi[twoI]   -= (dN_dx[count1][ii]*stre[0] + dN_dy[count1][ii]*stre[3]) ;
//           resi[twoIp1] -= (dN_dx[count1][ii]*stre[3] + dN_dy[count1][ii]*stre[1]) ;

           resi[twoI]   -= (Bbar[ii][0][0]*stre[0] + Bbar[ii][1][0]*stre[1] + Bbar[ii][2][0]*stre[2] + Bbar[ii][3][0]*stre[3]) ;
           resi[twoIp1] -= (Bbar[ii][0][1]*stre[0] + Bbar[ii][1][1]*stre[1] + Bbar[ii][2][1]*stre[2] + Bbar[ii][3][1]*stre[3]) ;

        }

        // 2.) Geometric Contribution
        //==============================================


        for(ii=0;ii<nlbf;ii++)
        {
           bc[0][0] = (Bbar[ii][0][0] * stre[0] + Bbar[ii][3][0] * stre[3] );
           bc[0][1] = (Bbar[ii][0][0] * stre[3] + Bbar[ii][3][0] * stre[1] );
           bc[0][2] =  Bbar[ii][1][0] * stre[3] ;
           bc[0][3] =  Bbar[ii][1][0] * stre[1] ;
           bc[0][4] =  Bbar[ii][2][0] * stre[2] ;

           bc[1][0] =  Bbar[ii][0][1] * stre[0] ;
           bc[1][1] =  Bbar[ii][0][1] * stre[3] ;
           bc[1][2] = (Bbar[ii][3][1] * stre[0] + Bbar[ii][1][1] * stre[3] );
           bc[1][3] = (Bbar[ii][3][1] * stre[3] + Bbar[ii][1][1] * stre[1] );
           bc[1][4] =  Bbar[ii][2][1] * stre[2] ;

           twoI   = 2*ii;
           twoIp1 = twoI+1;

           for(jj=0;jj<nlbf;jj++)
           {
              twoJ   = 2*jj;
              twoJp1 = twoJ+1;

              stiffness_local[twoI][twoJ]     +=  (bc[0][0] * Bbar[jj][0][0] + bc[0][1] * Bbar[jj][3][0] + bc[0][3] * Bbar[jj][1][0] + bc[0][4] * Bbar[jj][2][0]) ;
              stiffness_local[twoI][twoJp1]   +=  (bc[0][0] * Bbar[jj][0][1] + bc[0][2] * Bbar[jj][3][1] + bc[0][3] * Bbar[jj][1][1] + bc[0][4] * Bbar[jj][2][1]) ;
              stiffness_local[twoIp1][twoJ]   +=  (bc[1][0] * Bbar[jj][0][0] + bc[1][1] * Bbar[jj][3][0] + bc[1][3] * Bbar[jj][1][0] + bc[1][4] * Bbar[jj][2][0]) ;
              stiffness_local[twoIp1][twoJp1] +=  (bc[1][0] * Bbar[jj][0][1] + bc[1][2] * Bbar[jj][3][1] + bc[1][3] * Bbar[jj][1][1] + bc[1][4] * Bbar[jj][2][1]) ;
           }
        }




       // 3.) COMPUTE THE EXTRA TERMS
       //==============================================

//             fact = ( stre[0] + stre[1] + stre[2] )/9.0 ;

             fact = ( stre[0] + stre[1] )/9.0 ;

             fact9 = 3.0 * fact ;

             for(ii=0;ii<nlbf;ii++)
             {
                fact1 = ( stre[0] * dN_dx[count1][ii] + stre[3] * dN_dy[count1][ii] )/3.0 + fact * gmN[0][ii] ;
                fact2 = ( stre[3] * dN_dx[count1][ii] + stre[1] * dN_dy[count1][ii] )/3.0 + fact * gmN[1][ii] ;

                fact5 = fact * g[0][ii] ;
                fact6 = fact * g[1][ii] ;

                fact7 = fact9 * dN_dx[count1][ii] ;
                fact8 = fact9 * dN_dy[count1][ii] ;

                twoI   = 2*ii;
                twoIp1 = twoI+1;

                for(jj=0;jj<nlbf;jj++)
                {
                   twoJ   = 2*jj;
                   twoJp1 = twoJ+1;

                   fact3 = ( stre[0] * dN_dx[count1][jj] + stre[3] * dN_dy[count1][jj] )/3.0;
                   fact4 = ( stre[3] * dN_dx[count1][jj] + stre[1] * dN_dy[count1][jj] )/3.0;

                   stiffness_local[twoI][twoJ]     += (gmN[0][ii] * fact3 + fact1 * gmN[0][jj] - fact5 * g[0][jj] + fact7 * dN_dx[count1][jj]) ;
                   stiffness_local[twoI][twoJp1]   += (gmN[0][ii] * fact4 + fact1 * gmN[1][jj] - fact5 * g[1][jj] + fact8 * dN_dx[count1][jj]) ;
                   stiffness_local[twoIp1][twoJ]   += (gmN[1][ii] * fact3 + fact2 * gmN[0][jj] - fact6 * g[0][jj] + fact7 * dN_dy[count1][jj]) ;
                   stiffness_local[twoIp1][twoJp1] += (gmN[1][ii] * fact4 + fact2 * gmN[1][jj] - fact6 * g[1][jj] + fact8 * dN_dy[count1][jj]) ;

                }
             }


             for(ii=0;ii<nsize;ii++)
             {
                for(jj=0;jj<nsize;jj++)
                {
                   stiffness_local[ii][jj]  +=  fact * ( g2[ii][jj] - 3.0 * g3[ii][jj] );
                }
             }

       count++;
       count1++;
       ll += nivGP;
    }
  }
*/

  return 0;
}





int NurbsElem2DStructFbarSolid::calcStiffnessAndResidual2()
{
/*
   if(!finite || axsy)
   {
     cout << " Small strains OR AXISYMMETRIC cases are not supported by NurbsElem2DStructFbarSolid  " << endl;
     return -1;
   }

   int   err, isw, count, count1, ll, ii, jj, kk, mm, twoI, twoJ, twoIp1, twoJp1, index, cnt, gp1, gp2;

   int  p  = surf0->p,
        q  = surf0->q;

   double F[nGP][4], Flow[4], detF[nGP], detF0, F33=0.0, fact, fact1, fact2, fact3, fact4, stre[4], bc[2][4], alpha;

   double cc[4][4], Bbar[nlbf][3][2], g[2][nlbf], gmN[2][nlbf];

   double dvol, dvol0=0.0, utemp, vtemp, dt, J, volume[nGP], fact5, fact6, fact7, fact8, fact9;

   double N[nGP][nlbf], dN_dx[nGP][nlbf], dN_dy[nGP][nlbf];

   double udiff, vdiff, bb1, bb2, dummy;

   udiff = uvalues[1] - uvalues[0];
   vdiff = vvalues[1] - vvalues[0];

    // memory allocation for extra variables

    int  sizeM = p*q;

    double  Nbar[nGP][sizeM];


    MatrixXd  MAB(sizeM,sizeM), invMAB(sizeM,sizeM), Cmat1(sizeM,nlbf), Cmat2(sizeM,nlbf), Wmat1(sizeM,nlbf), Wmat2(sizeM,nlbf);

    double  detFbar, detJbar1[sizeM], detJbar2[sizeM], g2[nsize][nsize], g3[nsize][nsize];

    double  Cmat3[sizeM][nsize][nsize], Wmat3[sizeM][nsize][nsize], Cmat4[sizeM][nsize][nsize], Wmat4[sizeM][nsize][nsize]  ;


    MAB.setZero();
    invMAB.setZero();
    Cmat1.setZero();
    Cmat2.setZero();
    Wmat1.setZero();
    Wmat2.setZero();

   KNOTVECTOR  Ulow, Vlow;

   getLowerOrderKnotVector(surf1->U, surf1->p, Ulow);
   getLowerOrderKnotVector(surf1->V, surf1->q, Vlow);


    for(kk=0;kk<sizeM;kk++)
    {
       detJbar1[kk] = 0.0;
       detJbar2[kk] = 0.0;

       for(ii=0;ii<nsize;ii++)
       {
          for(jj=0;jj<nsize;jj++)
          {
             Cmat3[kk][ii][jj] = 0.0;
             Wmat3[kk][ii][jj] = 0.0;
             Cmat4[kk][ii][jj] = 0.0;
             Wmat4[kk][ii][jj] = 0.0;
          }
       }
    }

    for(ii=0;ii<4;ii++)
    {
       for(jj=0;jj<4;jj++)
       {
         cc[ii][jj] = 0.0;
       }
    }

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

           NurbsShapeFunctions2DAlg55(surf1, startindex[0], startindex[1], 1, dN_dx[count1], dN_dy[count1], F[count1], detF[count1]);

           NurbsShapeFunctions2DAlg11(surf1, startindex[0], startindex[1], knotsAtGPs[index], knotsAtGPs[index+1], N[count1], dN_dx[count1], dN_dy[count1], J);

           dvol0 = J * gaussweights[count1] * thick;

           volume[count1] = dvol0;


           dummy = pow(detF[count1], 0.5);

           // compute Mass matrix and (detF)^(1/2)

           for(ii=0;ii<sizeM;ii++)
           {
               fact = Nbar[count1][ii] * dvol0;

               detJbar1[ii] += fact * dummy;

               for(jj=0;jj<sizeM;jj++)
               {
                  MAB(ii,jj) += fact * Nbar[count1][jj];
               }
           }

           // compute C matrix

           for(ii=0;ii<sizeM;ii++)
           {
              fact = dvol0 * Nbar[count1][ii] * dummy;

              for(jj=0;jj<nlbf;jj++)
              {
                 Cmat1(ii,jj) += dN_dx[count1][jj] * fact;
                 Cmat2(ii,jj) += dN_dy[count1][jj] * fact;
              }
           }


           for(kk=0;kk<sizeM;kk++)
           {
              fact = dvol0 * Nbar[count1][kk] * dummy;

              for(ii=0;ii<nlbf;ii++)
              {
                 twoI   = 2*ii;
                 twoIp1 = twoI+1;

                 fact1 = fact * dN_dx[count1][ii];
                 fact2 = fact * dN_dy[count1][ii];

                 for(jj=0;jj<nlbf;jj++)
                 {
                    twoJ   = 2*jj;
                    twoJp1 = twoJ+1;

                    Cmat3[kk][twoI][twoJ]     +=  ( fact1 * dN_dx[count1][jj] );
                    Cmat3[kk][twoI][twoJp1]   +=  ( fact1 * dN_dy[count1][jj] );
                    Cmat3[kk][twoIp1][twoJ]   +=  ( fact2 * dN_dx[count1][jj] );
                    Cmat3[kk][twoIp1][twoJp1] +=  ( fact2 * dN_dy[count1][jj] );

                    Cmat4[kk][twoI][twoJ]     +=  ( fact1 * dN_dx[count1][jj] );
                    Cmat4[kk][twoI][twoJp1]   +=  ( fact2 * dN_dx[count1][jj] );
                    Cmat4[kk][twoIp1][twoJ]   +=  ( fact1 * dN_dy[count1][jj] );
                    Cmat4[kk][twoIp1][twoJp1] +=  ( fact2 * dN_dy[count1][jj] );

                 }
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
         {
            for(jj=0;jj<nsize;jj++)
            {
               for(kk=0;kk<sizeM;kk++)
               {
                  for(ll=0;ll<sizeM;ll++)
                  {
                     Wmat3[kk][ii][jj] += invMAB(kk,ll) * Cmat3[ll][ii][jj];
                     Wmat4[kk][ii][jj] += invMAB(kk,ll) * Cmat4[ll][ii][jj];
                  }
               }
            }
         }

         for(kk=0;kk<sizeM;kk++)
         {
            detJbar2[kk] = 0.0;
            for(ll=0;ll<sizeM;ll++)
               detJbar2[kk] += invMAB(kk,ll) * detJbar1[ll];
         }





    for(ii=0;ii<nsize;ii++)
      stiffness_local[ii].zero();

    resi.zero();

    count  = 1;
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

         detFbar = 0.0;
         for(ii=0;ii<sizeM;ii++)
            detFbar += detJbar2[ii] * Nbar[count1][ii];

         alpha = detFbar/pow(detF[count1], 0.5);


         for(ii=0;ii<nlbf;ii++)
         {
            g[0][ii] = 0.0;
            g[1][ii] = 0.0;
            for(jj=0;jj<sizeM;jj++)
            {
               g[0][ii] += Wmat1(jj,ii) * Nbar[count1][jj];
               g[1][ii] += Wmat2(jj,ii) * Nbar[count1][jj];
            }
            g[0][ii] /= detFbar;
            g[1][ii] /= detFbar;
         }


         for(ii=0;ii<nsize;ii++)
         {
            for(jj=0;jj<nsize;jj++)
            {
               g2[ii][jj] = 0.0;
               g3[ii][jj] = 0.0;

               for(kk=0;kk<sizeM;kk++)
               {
                   g2[ii][jj] +=  ( Nbar[count1][kk] * Wmat3[kk][ii][jj] );
                   g3[ii][jj] +=  ( Nbar[count1][kk] * Wmat4[kk][ii][jj] );
               }

               g2[ii][jj] /= detFbar;
               g3[ii][jj] /= detFbar;
            }
         }


         // calculate Fbar

         for(ii=0;ii<4;ii++)
             F[count1][ii] = alpha * F[count1][ii];

         F33 = alpha ;



        dt = mpapTime.dt;
        matlib2d_(matDat, F[count1], &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, (Element*) this);
        if(err !=0)           return 1;

//        if(elenum ==0 )        cout << '\t' << stre[0] << '\t' << stre[1] << '\t' << stre[3] << '\t' << stre[2] << endl; cout << endl;

        dvol0 = dvol0 * F33;


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
           gmN[0][ii] = g[0][ii] - dN_dx[count1][ii];
           gmN[1][ii] = g[1][ii] - dN_dy[count1][ii];

           bb1 = gmN[0][ii]/2.0;
           bb2 = gmN[1][ii]/2.0;

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

        // 1.) Material Contribution
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

           for(jj=0;jj<nlbf;jj++)
           {
              twoJ   = 2*jj;
              twoJp1 = twoJ+1;

              stiffness_local[twoI][twoJ]     +=  (bc[0][0] * Bbar[jj][0][0] + bc[0][1] * Bbar[jj][1][0] + bc[0][2] * Bbar[jj][2][0] ) ;
              stiffness_local[twoI][twoJp1]   +=  (bc[0][0] * Bbar[jj][0][1] + bc[0][1] * Bbar[jj][1][1] + bc[0][2] * Bbar[jj][2][1] ) ;
              stiffness_local[twoIp1][twoJ]   +=  (bc[1][0] * Bbar[jj][0][0] + bc[1][1] * Bbar[jj][1][0] + bc[1][2] * Bbar[jj][2][0] ) ;
              stiffness_local[twoIp1][twoJp1] +=  (bc[1][0] * Bbar[jj][0][1] + bc[1][1] * Bbar[jj][1][1] + bc[1][2] * Bbar[jj][2][1] ) ;

           }

//           resi[twoI]   -= (dN_dx[count1][ii]*stre[0] + dN_dy[count1][ii]*stre[3]) ;
//           resi[twoIp1] -= (dN_dx[count1][ii]*stre[3] + dN_dy[count1][ii]*stre[1]) ;

           resi[twoI]   -= (Bbar[ii][0][0]*stre[0] + Bbar[ii][1][0]*stre[1] + Bbar[ii][2][0]*stre[3] ) ;
           resi[twoIp1] -= (Bbar[ii][0][1]*stre[0] + Bbar[ii][1][1]*stre[1] + Bbar[ii][2][1]*stre[3] ) ;

        }

        // 2.) Geometric Contribution
        //==============================================


        for(ii=0;ii<nlbf;ii++)
        {
           bc[0][0] = (Bbar[ii][0][0] * stre[0] + Bbar[ii][2][0] * stre[3] );
           bc[0][1] = (Bbar[ii][0][0] * stre[3] + Bbar[ii][2][0] * stre[1] );
           bc[0][2] =  Bbar[ii][1][0] * stre[3] ;
           bc[0][3] =  Bbar[ii][1][0] * stre[1] ;

           bc[1][0] =  Bbar[ii][0][1] * stre[0] ;
           bc[1][1] =  Bbar[ii][0][1] * stre[3] ;
           bc[1][2] = (Bbar[ii][2][1] * stre[0] + Bbar[ii][1][1] * stre[3] );
           bc[1][3] = (Bbar[ii][2][1] * stre[3] + Bbar[ii][1][1] * stre[1] );

           twoI   = 2*ii;
           twoIp1 = twoI+1;

           for(jj=0;jj<nlbf;jj++)
           {
              twoJ   = 2*jj;
              twoJp1 = twoJ+1;

              stiffness_local[twoI][twoJ]     +=  (bc[0][0] * Bbar[jj][0][0] + bc[0][1] * Bbar[jj][2][0] + bc[0][3] * Bbar[jj][1][0] ) ;
              stiffness_local[twoI][twoJp1]   +=  (bc[0][0] * Bbar[jj][0][1] + bc[0][2] * Bbar[jj][2][1] + bc[0][3] * Bbar[jj][1][1] ) ;
              stiffness_local[twoIp1][twoJ]   +=  (bc[1][0] * Bbar[jj][0][0] + bc[1][1] * Bbar[jj][2][0] + bc[1][3] * Bbar[jj][1][0] ) ;
              stiffness_local[twoIp1][twoJp1] +=  (bc[1][0] * Bbar[jj][0][1] + bc[1][2] * Bbar[jj][2][1] + bc[1][3] * Bbar[jj][1][1] ) ;

           }
        }




       // 3.) COMPUTE THE EXTRA TERMS
       //==============================================

             fact = ( stre[0] + stre[1] )/4.0 ;

             fact9 = 2.0 * fact ;

             for(ii=0;ii<nlbf;ii++)
             {
                fact1 = ( stre[0] * dN_dx[count1][ii] + stre[3] * dN_dy[count1][ii] )/2.0 + fact * gmN[0][ii] ;
                fact2 = ( stre[3] * dN_dx[count1][ii] + stre[1] * dN_dy[count1][ii] )/2.0 + fact * gmN[1][ii] ;

                fact5 = fact * g[0][ii] ;
                fact6 = fact * g[1][ii] ;

                fact7 = fact9 * dN_dx[count1][ii] ;
                fact8 = fact9 * dN_dy[count1][ii] ;

                twoI   = 2*ii;
                twoIp1 = twoI+1;

                for(jj=0;jj<nlbf;jj++)
                {
                   twoJ   = 2*jj;
                   twoJp1 = twoJ+1;

                   fact3 = ( stre[0] * dN_dx[count1][jj] + stre[3] * dN_dy[count1][jj] )/2.0;
                   fact4 = ( stre[3] * dN_dx[count1][jj] + stre[1] * dN_dy[count1][jj] )/2.0;

                   stiffness_local[twoI][twoJ]     += (gmN[0][ii] * fact3 + fact1 * gmN[0][jj] - fact5 * g[0][jj] + fact7 * dN_dx[count1][jj]) ;
                   stiffness_local[twoI][twoJp1]   += (gmN[0][ii] * fact4 + fact1 * gmN[1][jj] - fact5 * g[1][jj] + fact8 * dN_dx[count1][jj]) ;
                   stiffness_local[twoIp1][twoJ]   += (gmN[1][ii] * fact3 + fact2 * gmN[0][jj] - fact6 * g[0][jj] + fact7 * dN_dy[count1][jj]) ;
                   stiffness_local[twoIp1][twoJp1] += (gmN[1][ii] * fact4 + fact2 * gmN[1][jj] - fact6 * g[1][jj] + fact8 * dN_dy[count1][jj]) ;

                }
             }


             for(ii=0;ii<nsize;ii++)
             {
                for(jj=0;jj<nsize;jj++)
                {
                   stiffness_local[ii][jj]  +=  fact * ( g2[ii][jj] - 2.0 * g3[ii][jj]);
                }
             }

       count++;
       count1++;
       ll += nivGP;
    }
  }
*/

  return 0;
}






int NurbsElem2DStructFbarSolid::calcStiffnessAndResidual3()
{
/*
   double F[4], detF=0.0, F33=0.0, fact, fact1, fact2, fact3, fact4,
          stre[4], cc[4][4], bc[2][4], cc1[4][4],
          F0[4], detF0, trc, cch[3];

   double rad0=0.0, rad=0.0, dvol=0.0, dvol0=0.0, r1dF33=0.0, r1drad0=0.0, utemp, vtemp, dt, J;

   int   err   = 0,
         isw   = 3,
         count = 1, count1 = 0, index,
         ll    = 0;

   int   gp1, gp2, ii, jj, kk, twoI, twoIp1, twoJ, twoJp1;

    int p = surf0->p, q = surf0->q;

    double N[nlbf], dN_dx[nlbf], dN_dy[nlbf];
    double N1[nlbf], dN_dx1[nlbf], dN_dy1[nlbf];

    Klocal.setZero();
    Flocal.setZero();


        // calculate F0 and detF0
        //--------------------------

        utemp = 0.5*(uvalues[0] + uvalues[1]);
        vtemp = 0.5*(vvalues[0] + vvalues[1]);

        NurbsShapeFunctions2DAlg11(surf0, startindex[0], startindex[1], utemp, vtemp, N1, dN_dx1, dN_dy1, J);
        NurbsShapeFunctions2DAlg55(surf1, startindex[0], startindex[1], 1, dN_dx1, dN_dy1, F, detF);

        if(finite)
            NurbsShapeFunctions2DAlg11(surf1, startindex[0], startindex[1], utemp, vtemp, N1, dN_dx1, dN_dy1, J);        

        trc = F[0] + F[3];
        detF0 = detF;


    double *gaussweights = &(surf0->gaussweights[0]);

   // loop over Gauss points
   for(gp2=0;gp2<nGP2;gp2++)
   {
      for(gp1=0;gp1<nGP1;gp1++)
      {
        index = count1*2;

        // COMPUTE SHAPE FUNCTIONS, THEIR DERIVATIVES IN UNDERFORMED CONFIGURATION

        NurbsShapeFunctions2DAlg11(surf0, startindex[0], startindex[1], knotsAtGPs[index], knotsAtGPs[index+1], N, dN_dx, dN_dy, J);

        // compute volume for current Gauss point in undeformed configuration
        // (for inertia and body forces)

        dvol0 = J * gaussweights[count1] * thick;

        // COMPUTE THE DEFORMATION GRADIENT 'F'
        NurbsShapeFunctions2DAlg55(surf1, startindex[0], startindex[1], 1, dN_dx, dN_dy, F, detF);

        // for finite deformation problems compute shape functions in deformed configuration
        if(finite)
        {
            NurbsShapeFunctions2DAlg11(surf1, startindex[0], startindex[1], knotsAtGPs[index], knotsAtGPs[index+1], N, dN_dx, dN_dy, J);        
        }

        dvol = J * gaussweights[count1] * thick;


        // ADJUST F33 fOR 2D PROBLEMS BASED ON THE ASSUMPTIONS OF PLANE STRESS/PLANE STRAIN/AXISYMMETRIC

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

        r1dF33 = 1.0/F33;


        // calculate Fbar

        if(finite)
        {
           fact = sqrt((detF0/detF));
           for(ii=0;ii<4;ii++)
             F[ii] = fact*F[ii];
        }
        else
        {
           fact = 0.5*(trc-F[0]-F[3]);
           F[0] += fact;
           F[3] += fact;
        }

        // COMPUTE MATERIAL RESPONSE

        dt = mpapTime.dt;
        matlib2d_(matDat, F, &F33, stre, cc1[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, (Element*) this);
        count++;
        count1++;
        ll += nivGP;


        for(ii=0;ii<4;ii++)
        {
           for(jj=0;jj<4;jj++)
           {
             cc[ii][jj] = cc1[ii][jj] ;
           }
        }

        cc[0][2] = cc1[0][3];
        cc[1][2] = cc1[1][3];

        cc[2][0] = cc1[3][0];
        cc[2][1] = cc1[3][1];

        cc[2][2] = cc1[3][3];
        cc[2][3] = cc1[3][3];
        cc[3][2] = cc1[3][3];


        if(err !=0)
          return 1;

        // MULTIPLY STRESS AND TANGENT TENSOR WITH VOLUME ELEMENT

        if(finite)
          dvol = dvol * F33;

        for(ii=0;ii<4;ii++)
        {
          stre[ii] *= dvol;
          for(jj=0;jj<4;jj++)
            cc[ii][jj] *= dvol;
        }



        //==============================================
        // CALCULATE TANGENT STIFFNESS
        //==============================================


       // correct tangent cc for Fbar forumation
            cch[0]  = 0.5 * (cc[0][0] + cc[0][1]);
            cch[1]  = 0.5 * (cc[1][0] + cc[1][1]);
            cch[2]  = 0.5 * (cc[2][0] + cc[2][1]);
            cc[0][0] = cc[0][0] - cch[0];
            cc[1][0] = cc[1][0] - cch[1];
            cc[2][0] = cc[2][0] - cch[2];
            cc[0][1] = cc[0][1] - cch[0];
            cc[1][1] = cc[1][1] - cch[1];
            cc[2][1] = cc[2][1] - cch[2];

       //   part 1. -- material part (not necessarily symmetric!!)  and   internal forces

       for(ii=0;ii<nlbf;ii++)
       {
           bc[0][0] = (dN_dx[ii] * cc[0][0] + dN_dy[ii] * cc[2][0]);
           bc[0][1] = (dN_dx[ii] * cc[0][1] + dN_dy[ii] * cc[2][1]);
           bc[0][2] = (dN_dx[ii] * cc[0][2] + dN_dy[ii] * cc[2][2]);
           bc[1][0] = (dN_dy[ii] * cc[1][0] + dN_dx[ii] * cc[2][0]);
           bc[1][1] = (dN_dy[ii] * cc[1][1] + dN_dx[ii] * cc[2][1]);
           bc[1][2] = (dN_dy[ii] * cc[1][2] + dN_dx[ii] * cc[2][2]);

           twoI = 2*ii;
           twoIp1 = twoI + 1;

           for(jj=0;jj<nlbf;jj++)
           {
              twoJ = 2*jj;
              twoJp1 = twoJ + 1;

              Klocal(twoI,twoJ)     +=  bc[0][0] * dN_dx[jj] + bc[0][2] * dN_dy[jj];
              Klocal(twoI,twoJp1)   +=  bc[0][1] * dN_dy[jj] + bc[0][2] * dN_dx[jj];
              Klocal(twoIp1,twoJ)   +=  bc[1][0] * dN_dx[jj] + bc[1][2] * dN_dy[jj];
              Klocal(twoIp1,twoJp1) +=  bc[1][1] * dN_dy[jj] + bc[1][2] * dN_dx[jj];
           }

           Flocal[twoI]   -= (dN_dx[ii]*stre[0] + dN_dy[ii]*stre[3]) ;
           Flocal[twoIp1] -= (dN_dx[ii]*stre[3] + dN_dy[ii]*stre[1]) ;
       }

        //   part 2. -- geometrical matrix  (if geometry nonlinear)

          if(finite)
          {
             for(ii=0;ii<nlbf;ii++)
             {
                fact1 = (dN_dx[ii] * stre[0] + dN_dy[ii] * stre[3]) ;
                fact2 = (dN_dx[ii] * stre[3] + dN_dy[ii] * stre[1]) ;

                twoI = 2*ii;
                twoIp1 = twoI + 1;

                for(jj=0;jj<nlbf;jj++)
                {
                   fact = fact1 * dN_dx[jj] + fact2 * dN_dy[jj];

                   twoJ = 2*jj;
                   twoJp1 = twoJ + 1;

                   Klocal(twoI,twoJ)     += fact ;
                   Klocal(twoIp1,twoJp1) += fact ;
                }
             }
          }

        //   derivative w.r.t shape functions of centroid

             for(ii=0;ii<nlbf;ii++)
             {
                fact1 = (dN_dx[ii] * cch[0] + dN_dy[ii] * cch[2]) ;
                fact2 = (dN_dx[ii] * cch[2] + dN_dy[ii] * cch[1]) ;

                twoI = 2*ii;
                twoIp1 = twoI + 1;

                for(jj=0;jj<nlbf;jj++)
                {
                    twoJ = 2*jj;
                    twoJp1 = twoJ + 1;

                    Klocal(twoI,twoJ)     +=  fact1 * dN_dx1[jj] ;
                    Klocal(twoI,twoJp1)   +=  fact1 * dN_dy1[jj] ;
                    Klocal(twoIp1,twoJ)   +=  fact2 * dN_dx1[jj] ;
                    Klocal(twoIp1,twoJp1) +=  fact2 * dN_dy1[jj] ;
                }
             }
    }
  }
*/

  return 0;
}



int NurbsElem2DStructFbarSolid::calcInternalForces()
{

  return 0;
}




int NurbsElem2DStructFbarSolid::calcOutput(double u1, double v1)
{

  return 0;
}







void NurbsElem2DStructFbarSolid::discreteContourplot(int vartype, int varindex, int index, int nCol, double umin, double umax)
{
  if(index > nivGP)
  {
     cout << '\t' << " Error in NurbsElem2DStructFbarSolid::contourplot " << endl;
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

  du = (uvalues[1] - uvalues[0])/nGP1;
  dv = (vvalues[1] - vvalues[0])/nGP2;

  ListArray<EPOINT> S1;
  S1.setDim( (nGP1+1)*(nGP2+1) );

  int count=0;
  vv = vvalues[0];

  if(finite)
  {
     for(int jj=0;jj<=nGP2;jj++)
     {
        uu = uvalues[0];
        for(int ii=0;ii<=nGP1;ii++)
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
     for(int jj=0;jj<=nGP2;jj++)
     {
        uu = uvalues[0];
        for(int ii=0;ii<=nGP1;ii++)
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

  int n1, ind1, ind2;
  n1 = 3;
  count=0;

  for(int jj=0;jj<nGP2;jj++)
  {
      ind1 = (nGP1+1)*jj;
      ind2 = (nGP1+1)*(jj+1);
      for(int ii=0;ii<nGP1;ii++)
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
          //plot.triangleContourPlot(x1, x2, x3, u1, u1, u1, umin, umax, nCol);

          // contour plot for 2nd triangle
          //plot.triangleContourPlot(x1, x3, x4, u1, u1, u1, umin, umax, nCol);

          count++;
      }
  }

  return;
}








void NurbsElem2DStructFbarSolid::projectToKnots(bool extrapolateFlag, int vartype, int varindex, int index)
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




void NurbsElem2DStructFbarSolid::projectStress(int varindex, double* outval)
{
/*
   int nivEL = nGP * nivGP;
   for(int ii=0;ii<nivEL;ii++)
     intVar2[ii] = intVar1[ii];

   double F[4], detF=0.0, F33=0.0, fact, fact1, fact2, fact3, fact4,
          stre[4], cc[4][4], bc[2][4], cch[3];

   double rad0=0.0, rad=0.0, dvol=0.0, dvol0=0.0, r1dF33=0.0, r1drad0=0.0, utemp, vtemp, dt, J;

   int   err   = 0,
         isw   = 3,
         count = 1, count1 = 0, index,
         ll    = 0;

    int p = surf0->p, q = surf0->q, cnt=0;

    double N[nlbf], dN_dx[nlbf], dN_dy[nlbf];

    int ntmp1 = (nGP1-1), ntmp2 = (nGP2-1), ntmp = ntmp1*ntmp2;

    double udiff, vdiff, Jbar, Jbar1, detF0[ntmp], NNtmp[ntmp];

    double N1[nlbf], dN_dx1[ntmp][nlbf], dN_dy1[ntmp][nlbf], dN_dx2[nlbf], dN_dy2[nlbf];


    vector<double> gpts1, gpts2, gwts1, gwts2;

    // get the Gauss Points and Weights in the first direction

    getGaussPoints1D(ntmp1, gpts1, gwts1);
    getGaussPoints1D(ntmp2, gpts2, gwts2);

    udiff = uvalues[1] - uvalues[0];
    vdiff = vvalues[1] - vvalues[0];

    // calculate F, detF and gradient matrix at sampling points
    //--------------------------

    // loop over sampling points
    for(int gp2=0;gp2<ntmp2;gp2++)
    {
       for(int gp1=0;gp1<ntmp1;gp1++)
       {
	        utemp = uvalues[0] + 0.5*udiff*(1 + gpts1[gp1]);
	        vtemp = vvalues[0] + 0.5*vdiff*(1 + gpts2[gp2]);

	        NurbsShapeFunctions2DAlg11(surf0, startindex[0], startindex[1], utemp, vtemp, N1, dN_dx1[cnt], dN_dy1[cnt], J);
	        NurbsShapeFunctions2DAlg55(surf1, startindex[0], startindex[1], 1, dN_dx1[cnt], dN_dy1[cnt], F, detF);

	        detF0[cnt] = detF;

	        cnt++;
	   }
    }



    double *gaussweights = &(surf0->gaussweights[0]);

    for(int ii=0;ii<nsize;ii++)
      stiffness_local[ii].zero();

    resi.zero();


   // loop over Gauss points
   for(int gp2=0;gp2<nGP2;gp2++)
   {
      for(int gp1=0;gp1<nGP1;gp1++)
      {
        index = count1*2;

        // COMPUTE SHAPE FUNCTIONS, THEIR DERIVATIVES IN UNDERFORMED CONFIGURATION

        NurbsShapeFunctions2DAlg11(surf0, startindex[0], startindex[1], knotsAtGPs[index], knotsAtGPs[index+1], N, dN_dx, dN_dy, J);

        // compute volume for current Gauss point in undeformed configuration
        // (for inertia and body forces)

        dvol0 = J * gaussweights[count1] * thick;

        //   for axisymmetric problems compute radius
        if(axsy)
        {
           rad  = rad0 = 0.0 ;

           int loc_num=0;
           for(int jj=0; jj<=surf0->q; jj++)
           {
              for(int ii=0; ii<=surf0->p; ii++)
              {
                 rad0 += (N[loc_num] * surf0->Pw[startindex[0]+ii][startindex[1]+jj].CalcEuclid().x);
                 rad  += (N[loc_num] * surf1->Pw[startindex[0]+ii][startindex[1]+jj].CalcEuclid().x);
                 loc_num++;
              }
           }
           dvol0 *= (twoPI * rad0);
           r1drad0 = 1.0 / rad0;
        }

        // COMPUTE THE DEFORMATION GRADIENT 'F'
        NurbsShapeFunctions2DAlg55(surf1, startindex[0], startindex[1], 1, dN_dx, dN_dy, F, detF);

        // for finite deformation problems compute shape functions in deformed configuration
        if(finite)
        {
            NurbsShapeFunctions2DAlg11(surf1, startindex[0], startindex[1], knotsAtGPs[index], knotsAtGPs[index+1], N, dN_dx, dN_dy, J);
        }

        // ADJUST F33 fOR 2D PROBLEMS BASED ON THE ASSUMPTIONS OF PLANE STRESS/PLANE STRAIN/AXISYMMETRIC

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

        r1dF33 = 1.0/F33;


        // calculate Fbar

        if(finite)
        {
		double NNtmp1[ntmp1], NNtmp2[ntmp2];
		Bernstein_Basis(ntmp1-1, knotsAtGPs[index]-uvalues[0], NNtmp1);
		Bernstein_Basis(ntmp2-1, knotsAtGPs[index+1]-vvalues[0], NNtmp2);

		cnt = 0;
		for(int a2=0;a2<ntmp2;a2++)
		{
			for(int a1=0;a1<ntmp1;a1++)
			{
				NNtmp[cnt++] = NNtmp1[a1] * NNtmp2[a2];
			}
		}

		Jbar1 = 0.0;

		for(int kk=0;kk<ntmp;kk++)
		    Jbar1 += detF0[kk] * NNtmp[kk];

            fact = sqrt((Jbar1/detF));
            for(int ii=0;ii<4;ii++)
              F[ii] = fact*F[ii];
        }

        // COMPUTE MATERIAL RESPONSE

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
           cout << '\t' << "    NurbsElem2DStructFbarSolid::projectStress .... : Error in 'varindex' " << endl;
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




void NurbsElem2DStructFbarSolid::projectStrain(int vartype, int varindex, double* outval)
{
/*
   double F[4], detF=0.0, F33=0.0, stre[4], cc[4][4];

   double rad0=0.0, rad=0.0, r1dF33=0.0, r1drad0=0.0, J=0.0, dt;

   int   err    = 0,
         isw    = 3,
         count  = 1, 
         count1 = 0, index, ll = 0, ii, jj;

   double *gaussweights = &(surf0->gaussweights[0]);

   int p = surf0->p, q = surf0->q;

   double N[nlbf], dN_dx[nlbf], dN_dy[nlbf];

   int nivEL = nGP * nivGP;
   for(ii=0;ii<nivEL;ii++)
     intVar2[ii] = intVar1[ii];

   bool intVarFlag = false;
   if( nivGP > 0)
      intVarFlag = true;



   // loop over Gauss points
   for(int gp2=0;gp2<nGP2;gp2++)
   {
      for(int gp1=0;gp1<nGP1;gp1++)
      {
        // COMPUTE SHAPE FUNCTIONS, THEIR DERIVATIVES IN UNDERFORMED CONFIGURATION

        index = count1*2;

        NurbsShapeFunctions2DAlg11(surf0, startindex[0], startindex[1], knotsAtGPs[index], knotsAtGPs[index+1], N, dN_dx, dN_dy, J);

        if(J < 0.0)
        {
           cerr << "        NurbsElem2DStructFbarSolid::projectStress.......Negative Jacobian   " << endl;
           return ;
        }

        //   for axisymmetric problems compute radius
        if(axsy)
        {
           rad  = rad0 = 0.0 ;

           int loc_num=0, temp1;
           for(int jj=0; jj<=q; jj++)
           {
              temp1 = startindex[1]+jj;
              for(int ii=0; ii<=p; ii++)
              {
                 rad0 += (N[loc_num] * surf0->Pw[startindex[0]+ii][temp1].CalcEuclid().x);
                 rad  += (N[loc_num] * surf1->Pw[startindex[0]+ii][temp1].CalcEuclid().x);
                 loc_num++;
              }
           }
        }

        // COMPUTE THE DEFORMATION GRADIENT 'F'
        NurbsShapeFunctions2DAlg55(surf1, startindex[0], startindex[1], 1, dN_dx, dN_dy, F, detF);

        // for finite deformation problems compute shape functions in deformed configuration
        if(finite)
        {
            NurbsShapeFunctions2DAlg11(surf1, startindex[0], startindex[1], knotsAtGPs[index], knotsAtGPs[index+1], N, dN_dx, dN_dy, J);
        }

        // ADJUST F33 fOR 2D PROBLEMS BASED ON THE ASSUMPTIONS OF PLANE STRESS/PLANE STRAIN/AXISYMMETRIC

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

        // COMPUTE MATERIAL RESPONSE

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










void NurbsElem2DStructFbarSolid::projectIntVar(int index, double* outval)
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





void  NurbsElem2DStructFbarSolid::AssembleElementMatrix(int index, SparseMatrixXd& mtx, int start1, int start2)
{
    int  ii, jj, nn=0, aa, bb, ind, *tt1;

    tt1 = &(surf0->LM[elenum][0]);

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



