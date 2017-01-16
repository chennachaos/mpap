#include <iostream>
#include <math.h>
#include "NurbsSURFACE.h"
#include "Debug.h"
#include "NurbsShapeFunctions.h"

using namespace std;



void NurbsShapeFunctions2DAlg2(NurbsSURFACE* surf1, int ni, int nj, double u_tilde, double v_tilde, double* NN, double& J, double* dircos)
{
   short p = surf1->p, q = surf1->q, ii, jj, loc_num = 0;

   int ind1 = ni+p, ind2 = ind1+1, ind3 = nj+q, ind4 = ind3+1;

   double u, v, dx_du, dx = 0.0, dy = 0.0 ;

   int ROWS = 2;

   double** ders1 = new double*[ROWS];

  EPOINT *EP;

  if(CompareDoubles(u_tilde,-1.0))  // side #4
  {
     v = ((surf1->V[ind4]-surf1->V[ind3])*v_tilde + (surf1->V[ind4]+surf1->V[ind3])) / 2.0;

     for(ii=0;ii<ROWS;ii++)
       ders1[ii] = new double[q+1];

     DersBasisFuns(&(surf1->V[0]), surf1->V.n, q, v, 1, ders1);

    vector<double>  dR_du(q+1);

    double sum_tot=0.0, sum1=0.0, wt;

    for(ii=0; ii<=q; ii++)
    {
      wt = surf1->Pw[ni][nj+ii].w ;

      NN[ii] = ders1[0][ii] * wt;
      
      sum_tot += NN[ii] ;

      dR_du[ii] = ders1[1][ii] * wt;
      
      sum1 += dR_du[ii];
    }

    double fact = sum_tot*sum_tot;

    for(loc_num = 0;loc_num<=q; loc_num++)
    {
      dR_du[loc_num] = (dR_du[loc_num]*sum_tot - NN[loc_num]*sum1)/fact;

      NN[loc_num] = NN[loc_num]/sum_tot;
    }

    dx = 0.0; dy = 0.0 ;
    for(jj=0; jj<=q; jj++)
    {
      EP = &(surf1->PP[ni][nj+jj]);

      dx += (EP->x * dR_du[jj]) ;
      dy += (EP->y * dR_du[jj]) ;
    }

     dx_du = sqrt(dx*dx + dy*dy);

     dircos[0] = dy/dx_du;
     dircos[1] = -dx/dx_du;

     J = dx_du * (surf1->V[ind4]-surf1->V[ind3])/2.0;
  }

  if(CompareDoubles(u_tilde,1.0))  // side #2
  {
     v = ((surf1->V[ind4]-surf1->V[ind3])*v_tilde + (surf1->V[ind4]+surf1->V[ind3])) / 2.0;

     for(ii=0;ii<ROWS;ii++)
       ders1[ii] = new double[q+1];

     DersBasisFuns(&(surf1->V[0]), surf1->V.n, q, v, 1, ders1);

    vector<double>  dR_du(q+1);

    double sum_tot=0.0, sum1=0.0, wt;

    for(ii=0; ii<=q; ii++)
    {
      wt = surf1->Pw[ind1][nj+ii].w ;

      NN[ii] = ders1[0][ii] * wt;
      
      sum_tot += NN[ii] ;

      dR_du[ii] = ders1[1][ii] * wt;
      
      sum1 += dR_du[ii];
    }

    double fact = sum_tot*sum_tot;

    for(loc_num = 0;loc_num<=q; loc_num++)
    {
      dR_du[loc_num] = (dR_du[loc_num]*sum_tot - NN[loc_num]*sum1)/fact;

      NN[loc_num] = NN[loc_num]/sum_tot;
    }

    dx = 0.0; dy = 0.0 ;
    for(jj=0; jj<=q; jj++)
    {
      EP = &(surf1->PP[ind1][nj+jj]);

      dx += (EP->x * dR_du[jj]);
      dy += (EP->y * dR_du[jj]);
    }

    dx_du = sqrt(dx*dx + dy*dy);

    dircos[0] = dy/dx_du;
    dircos[1] = -dx/dx_du;

    J = dx_du * (surf1->V[ind4]-surf1->V[ind3])/2.0;
  }

  if(CompareDoubles(v_tilde,-1.0))  // side #1
  {
     u = ((surf1->U[ind2]-surf1->U[ind1])*u_tilde + (surf1->U[ind2]+surf1->U[ind1])) / 2.0;

     for(ii=0;ii<ROWS;ii++)
       ders1[ii] = new double[p+1];

    DersBasisFuns(&(surf1->U[0]), surf1->U.n, p, u, 1, ders1);
    
    vector<double>  dR_du(p+1);

    double sum_tot=0.0, sum1=0.0, wt;

    for(ii=0; ii<=p; ii++)
    {
      wt = surf1->Pw[ni+ii][nj].w ;
      
      NN[ii] = ders1[0][ii] * wt;
      
      sum_tot += NN[ii] ;

      dR_du[ii] = ders1[1][ii] * wt;
      
      sum1 += dR_du[ii];
    }

    double fact = sum_tot*sum_tot;

    for(loc_num = 0;loc_num<=p; loc_num++)
    {
      dR_du[loc_num] = (dR_du[loc_num]*sum_tot - NN[loc_num]*sum1)/fact;

      NN[loc_num] = NN[loc_num]/sum_tot;
    }

    dx = 0.0; dy = 0.0 ;
    for(ii=0; ii<=p; ii++)
    {
      EP = &(surf1->PP[ni+ii][nj]);

      dx += ( EP->x * dR_du[ii]) ;
      dy += ( EP->y * dR_du[ii]) ;
    }

    dx_du = sqrt(dx*dx + dy*dy);

    dircos[0] = dy/dx_du;
    dircos[1] = -dx/dx_du;

    J = dx_du * (surf1->U[ind2]-surf1->U[ind1])/2.0;
  }

  if(CompareDoubles(v_tilde,1.0))  // side #3
  {
    u = ((surf1->U[ind2]-surf1->U[ind1])*u_tilde + (surf1->U[ind2]+surf1->U[ind1])) / 2.0;

    for(ii=0;ii<ROWS;ii++)
      ders1[ii] = new double[p+1];

    DersBasisFuns(&(surf1->U[0]), surf1->U.n, p, u, 1, ders1);

    vector<double>  dR_du(p+1);

    double sum_tot=0.0, sum1=0.0, wt;

    for(ii=0; ii<=p; ii++)
    {
      wt = surf1->Pw[ni+ii][ind3].w ;
      
      NN[ii] = ders1[0][ii] * wt;
      
      sum_tot += NN[ii] ;

      dR_du[ii] = ders1[1][ii] * wt;
      
      sum1 += dR_du[ii];
    }

    double fact = sum_tot*sum_tot;

    for(loc_num = 0;loc_num<=p; loc_num++)
    {
      dR_du[loc_num] = (dR_du[loc_num]*sum_tot - NN[loc_num]*sum1)/fact;

      NN[loc_num] = NN[loc_num]/sum_tot;
    }

    dx = 0.0; dy = 0.0 ;
    for(ii=0; ii<=p; ii++)
    {
      EP = &(surf1->PP[ni+ii][ind3]);

      dx += ( EP->x * dR_du[ii]) ;
      dy += ( EP->y * dR_du[ii]) ;
    }


    /*
    dx = 0.0; dy = 0.0 ;
    for(ii=0; ii<=p; ii++)
    {
      NN[ii] = ders1[0][ii] ;
      
      EP = &(surf1->PP[ni+ii][ind3]);

      dx += ( EP->x * ders1[1][ii]) ;
      dy += ( EP->y * ders1[1][ii]) ;
    }
    */

    dx_du = sqrt(dx*dx + dy*dy);

    dircos[0] = dy/dx_du;
    dircos[1] = -dx/dx_du;

    J = dx_du * (surf1->U[ind2]-surf1->U[ind1])/2.0;
  }

  for(ii=0;ii<ROWS;ii++)
    delete [] ders1[ii];

  delete [] ders1;

  return;
}
