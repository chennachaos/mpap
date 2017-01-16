#include <iostream>
#include <math.h>
#include "Debug.h"
#include "NurbsShapeFunctions.h"

using namespace std;


void NurbsShapeFunctions2DAlg55(NurbsSURFACE* surf1, int ni, int nj, int confflag, double* dN_dx, double* dN_dy, double* F, double& detF)
{
   // CALCULATION OF THE DEFORMATION GRADIENT F

      /////////////////////////////////////////////
      //                                         //
      //    F(1,1) = F[0]   //   F(1,2) = F[1]   //
      //                                         //
      //    F(2,1) = F[2]   //   F(2,2) = F[3]   //
      //                                         //
      /////////////////////////////////////////////

      //   confflag = 1 --> coordiantes from current configuration
      //                    shape function derivatives w.r.t reference configuration

      //   surf1 is updated NurbsSURFACE ( SurfaceResult[] in IsogeometricFEM )
      //   shpfnc1 is calculated using original NurbsSURFACE (SurfaceListFinal[] in IsogeometricFEM ) 
      //

      //   confflag = 2 ==> coordiantes from reference configuration
      //                    shape function derivatives w.r.t current configuration

      //   surf1 is original NurbsSURFACE (SurfaceListFinal[] in IsogeometricFEM )
      //   shpfnc1 is calculated using updated NurbsSURFACE ( SurfaceResult[] in IsogeometricFEM )
      //

      short ii, jj, loc_num;

      EPOINT *EP;

      F[0] = F[1] = F[2] = F[3] = 0.0;

      int temp;
      if( confflag == 1) // calculation of F
      {
         loc_num = 0;
         for(jj=0; jj<=surf1->q; jj++)
         {
            temp = nj+jj;
            for(ii=0; ii<=surf1->p; ii++)
            {
               //EP = surf1->Pw[ni+ii][temp].CalcEuclid();
               EP = &( surf1->PP[ni+ii][temp] );

               F[0] += EP->x * dN_dx[loc_num];
               F[2] += EP->x * dN_dy[loc_num];
               F[1] += EP->y * dN_dx[loc_num];
               F[3] += EP->y * dN_dy[loc_num];

               loc_num++;
            }
         }

         detF = F[0]*F[3] - F[1]*F[2];
      }
      else //confflag == 2.  calculation of F^-1
      {
         double Finv[4] = {0.0, 0.0, 0.0, 0.0};

         loc_num = 0;
         for(jj=0; jj<=surf1->q; jj++)
         {
            for(ii=0; ii<=surf1->p; ii++)
            {
               EP = &(surf1->PP[ni+ii][nj+jj] );

               Finv[0] += EP->x * dN_dx[loc_num];
               Finv[2] += EP->x * dN_dy[loc_num];
               Finv[1] += EP->y * dN_dx[loc_num];
               Finv[3] += EP->y * dN_dy[loc_num];
               loc_num++;
            }
         }

         detF = 1.0 / (Finv[0]*Finv[3] - Finv[1]*Finv[2]);

         F[0] =   Finv[3] * detF ;
         F[1] = - Finv[1] * detF ;
         F[2] = - Finv[2] * detF ;
         F[3] =   Finv[0] * detF ;
      }


  return;
}
