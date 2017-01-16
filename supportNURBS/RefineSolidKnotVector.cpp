#include <iostream>
#include <math.h>
#include "NurbsUtilitiesSOLID.h"

using namespace std;


void RefineSolidKnotVector(NurbsSOLID* solid1, int repamID, int* tt, NurbsSOLID* solid2)
{
   int   t1, t2, t3, ii, jj, kk;
   t1 = tt[0];
   t2 = tt[1];
   t3 = tt[2];

   solid2->Pw = solid1->Pw;
   solid2->U = solid1->U;
   solid2->V = solid1->V;
   solid2->W = solid1->W;
   solid2->p = solid1->p;
   solid2->q = solid1->q;
   solid2->r = solid1->r;
   
//   cout << t1 << '\t' << t2 << '\t' << t3 << endl;
   
   if(t1 == 0 && t2 == 0 && t3 == 0)
     return;


   KNOTVECTOR  U, V, W;
   DEGREE  p, q, r;

   U = solid1->U;
   V = solid1->V;
   W = solid1->W;
   p = solid1->p;
   q = solid1->q;
   r = solid1->r;

      int  ngbf1, ngbf2, ngbf3, ngbf1new, ngbf2new, ngbf3new;
      bool  flag1, flag2, flag3;
      flag1 = flag2 = flag3 = false;

      ngbf1 = U.n - p - 1;
      ngbf2 = V.n - q - 1;
      ngbf3 = W.n - r - 1;
      
      ngbf1new = ngbf1;
      ngbf2new = ngbf2;
      ngbf3new = ngbf3;

//      cout << " old ngbf ... : " << ngbf1 << '\t' << ngbf2 << '\t' << ngbf3 << endl;

      KNOTVECTOR  XU, XV, XW;
      
          if(repamID == 1 || repamID == 4)
          {
             if(t1 > 1)
             {
                KNOTVECTOR  XU1;

                create_vector2(U, t1, XU1);
             
                XU.setDim(XU1.n-2);
                for(ii=0;ii<XU.n;ii++)
                  XU[ii] = XU1[ii+1];
             
                ngbf1new = U.n + XU.n - p - 1;
                
                flag1 = true;
             }
          }
          else
          {
             if(t1 > 0)
             {
                GenKnotVecForRefining(U, t1, XU);
                ngbf1new = U.n + XU.n - p - 1;
                flag1 = true;
             }
          }

          if(repamID == 2 || repamID == 4)
          {
              if(t2 > 1)
              {
                 KNOTVECTOR  XV1;

                 create_vector2(V, t2, XV1);

                 XV.setDim(XV1.n-2);
                 for(ii=0;ii<XV.n;ii++)
                   XV[ii] = XV1[ii+1];

                 ngbf2new = V.n + XV.n - q - 1;
                 flag2 = true;
              }
          }
          else
          {
             if(t2 > 0)
             {
                GenKnotVecForRefining(V, t2, XV);
                ngbf2new = V.n + XV.n - q - 1;
                flag2 = true;
             }
          }

          if(repamID == 3 || repamID == 4)
          {
              if(t3 > 1)
              {
                 KNOTVECTOR  XW1;

                 create_vector2(W, t3, XW1);

                 XW.setDim(XW1.n-2);
                 for(ii=0;ii<XW.n;ii++)
                   XW[ii] = XW1[ii+1];

                 ngbf3new = W.n + XW.n - r - 1;
                 flag3 = true;
              }
          }
          else
          {
              if(t3 > 0)
              {
                 GenKnotVecForRefining(W, t3, XW);
                 ngbf3new = W.n + XW.n - r - 1;
                 flag3 = true;
              }
          }

      //cout << XU << endl;
      //cout << XV << endl;
      //cout << XW << endl;

      //cout << " new ngbf ... : " << ngbf1new << '\t' << ngbf2new << '\t' << ngbf3new << endl;
      
      //cout << " flags ... : " << flag1 << '\t' << flag2 << '\t' << flag3 << endl;


       if( flag3 )
       {
           solid2->Pw.setDim(ngbf3new);
           for(kk=0;kk<ngbf3new;kk++)
           {
              solid2->Pw[kk].setDim(ngbf1);
              for(ii=0;ii<ngbf1;ii++)
                 solid2->Pw[kk][ii].setDim(ngbf2);
           }
           for(kk=0;kk<ngbf3;kk++)
           {
              for(ii=0;ii<ngbf1;ii++)
              {
                 for(jj=0;jj<ngbf2;jj++)
                   solid2->Pw[kk][ii][jj] = solid1->Pw[kk][ii][jj];
              }
           }

           NurbsSURFACE  surf_temp1(solid1->Pw[0], W, V, r, q);
           
           NurbsSURFACE  surf_temp2(solid1->Pw[0], W, V, r, q);
           
           CNET  Pw1;
           
           Pw1.setDim(ngbf3);
           for(kk=0;kk<ngbf3;kk++)
             Pw1[kk].setDim(ngbf2);

           
           for(ii=0;ii<ngbf1;ii++)
           {
               for(kk=0;kk<ngbf3;kk++)
               {
                  for(jj=0;jj<ngbf2;jj++)
                  {
                     Pw1[kk][jj] = solid2->Pw[kk][ii][jj];
                  }
               }

               surf_temp1.Pw = Pw1;
      
               RefineSurfKnotVector1D(&surf_temp1, XW, 1, &surf_temp2);

               for(kk=0;kk<ngbf3new;kk++)
               {
                  for(jj=0;jj<ngbf2;jj++)
                  {
                     solid2->Pw[kk][ii][jj] = surf_temp2.Pw[kk][jj];
                  }
               }
           }

           solid2->W = surf_temp2.U;
       }


       if( flag1 )
       {
           NurbsSURFACE  surf_temp1(solid2->Pw[0], U, V, p, q);
           
           NurbsSURFACE  surf_temp2(solid2->Pw[0], U, V, p, q);
           
           NurbsSURFACE  surf_temp3(solid2->Pw[0], U, V, p, q);

           for(kk=0;kk<ngbf3new;kk++)
           {
              surf_temp1.Pw = solid2->Pw[kk];
      
              RefineSurfKnotVector1D(&surf_temp1, XU, 1, &surf_temp2);

              solid2->Pw[kk] = surf_temp2.Pw;
           }

           solid2->U = surf_temp2.U;

           if(flag2)
           {
               for(kk=0;kk<ngbf3new;kk++)
               {
                   surf_temp2.Pw = solid2->Pw[kk];
      
                   RefineSurfKnotVector1D(&surf_temp2, XV, 2, &surf_temp3);

                   solid2->Pw[kk] = surf_temp3.Pw;
               }

               solid2->V = surf_temp3.V;
           }
       }
       else
       {
           if(flag2)
           {
               NurbsSURFACE  surf_temp1(solid2->Pw[0], U, V, p, q);
               
               NurbsSURFACE  surf_temp2(solid2->Pw[0], U, V, p, q);

               for(kk=0;kk<ngbf3new;kk++)
               {
                   surf_temp1.Pw = solid2->Pw[kk];
      
                   RefineSurfKnotVector1D(&surf_temp1, XV, 2, &surf_temp2);

                   solid2->Pw[kk] = surf_temp2.Pw;
               }

               solid2->V = surf_temp2.V;
           }
       }

  return;
}
