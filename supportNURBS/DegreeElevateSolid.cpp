#include <iostream>
#include <math.h>
#include "NurbsSOLID.h"

using namespace std;


void  DegreeElevateSolid(NurbsSOLID* solid1, int* tt, NurbsSOLID* solid2)
{
   int   t1, t2, t3, ii, jj, kk;

   solid2->Pw = solid1->Pw;
   solid2->U = solid1->U;
   solid2->V = solid1->V;
   solid2->W = solid1->W;
   solid2->p = solid1->p;
   solid2->q = solid1->q;
   solid2->r = solid1->r;
   

   DEGREE  p, q, r, ph, qh, rh;

   p = solid1->p;
   q = solid1->q;
   r = solid1->r;

   t1 = tt[0] - p;
   t2 = tt[1] - q;
   t3 = tt[2] - r;

   //cout << t1 << '\t' << t2 << '\t' << t3 << endl;
   
   if(t1 == 0 && t2 == 0 && t3 == 0)
     return;


   KNOTVECTOR  U, V, W;
   U = solid1->U;
   V = solid1->V;
   W = solid1->W;

   ph = p + t1;
   qh = q + t2;
   rh = r + t3;

      int  ngbf1, ngbf2, ngbf3, ngbf1new, ngbf2new, ngbf3new;
           
      ngbf1 = U.n - p - 1;
      ngbf2 = V.n - q - 1;
      ngbf3 = W.n - r - 1;
      
      ngbf1new = ngbf1;
      ngbf2new = ngbf2;
      ngbf3new = ngbf3;

      //cout << " old ngbf ... : " << ngbf1 << '\t' << ngbf2 << '\t' << ngbf3 << endl;

      KNOTVECTOR  XU, XV, XW;
      
      if(t1 > 0)
      {
         findunique(U, XU);
         ngbf1new = U.n + t1*XU.n - ph - 1;
      }
      if(t2 > 0)
      {
         findunique(V, XV);
         ngbf2new = V.n + t2*XV.n - qh - 1;
      }
      if(t3 > 0)
      {
         findunique(W, XW);
         ngbf3new = W.n + t3*XW.n - rh - 1;
      }           

      //cout << " new ngbf ... : " << ngbf1new << '\t' << ngbf2new << '\t' << ngbf3new << endl;


       if(t3 > 0)
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
      
               DegreeElevateSurf1D(&surf_temp1, t3, 1, &surf_temp2);

               for(kk=0;kk<ngbf3new;kk++)
               {
                  for(jj=0;jj<ngbf2;jj++)
                  {
                     solid2->Pw[kk][ii][jj] = surf_temp2.Pw[kk][jj];
                  }
               }
           }

           solid2->W = surf_temp2.U;
           solid2->r = rh;
       }


       if(t1 > 0)
       {
           NurbsSURFACE  surf_temp1(solid2->Pw[0], U, V, p, q);
           
           NurbsSURFACE  surf_temp2(solid2->Pw[0], U, V, p, q);
           
           NurbsSURFACE  surf_temp3(solid2->Pw[0], U, V, p, q);

           for(kk=0;kk<ngbf3new;kk++)
           {
              surf_temp1.Pw = solid2->Pw[kk];
      
              DegreeElevateSurf1D(&surf_temp1, t1, 1, &surf_temp2);

              solid2->Pw[kk] = surf_temp2.Pw;
           }

           solid2->U = surf_temp2.U;
           solid2->p = ph;

           if(t2 > 0)
           {
               for(kk=0;kk<ngbf3new;kk++)
               {
                   surf_temp2.Pw = solid2->Pw[kk];
      
                   DegreeElevateSurf1D(&surf_temp2, t2, 2, &surf_temp3);

                   solid2->Pw[kk] = surf_temp3.Pw;
               }

               solid2->V = surf_temp3.V;
               solid2->q = qh;
           }
       }
       else if(t1 == 0)
       {
           if(t2 > 0)
           {
               NurbsSURFACE  surf_temp1(solid2->Pw[0], U, V, p, q);
               
               NurbsSURFACE  surf_temp2(solid2->Pw[0], U, V, p, q);

               for(kk=0;kk<ngbf3new;kk++)
               {
                   surf_temp1.Pw = solid2->Pw[kk];
      
                   DegreeElevateSurf1D(&surf_temp1, t2, 2, &surf_temp2);

                   solid2->Pw[kk] = surf_temp2.Pw;
               }

               solid2->V = surf_temp2.V;
               solid2->q = qh;
           }
       }

  return;
}
