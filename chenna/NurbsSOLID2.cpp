/*=============================================================================
        File: NurbsSOLID.cpp
  Created by: Chennakesava Kadapa          (08 Jan 2011)
 Purpose    : Implementation file for the definitions of NURBS SOLID Class

 ============================================================================*/

#include <Eigen/Dense>

#include "NurbsSOLID.h"
#include "DataBlockTemplate.h"
#include <iomanip>

using namespace std;
using namespace Eigen;


extern MpapTime mpapTime;






void NurbsSOLID::ShapeFunctions(double u, double v, double w, double* NN)
{
   vector<double>  Nu(p+1), Nv(q+1), Nw(r+1);
   double  temp;

   BasisFuns(&(U[0]), U.n, p, u, &Nu[0]);
   BasisFuns(&(V[0]), V.n, q, v, &Nv[0]);
   BasisFuns(&(W[0]), W.n, r, w, &Nw[0]);

   int loc_num = 0, ii, jj, kk;

   for(kk=0;kk<=r;kk++)
   {
      for(jj=0;jj<=q;jj++)
      {
         temp = Nw[kk] * Nv[jj];
         for(ii=0;ii<=p;ii++)
         {
            NN[loc_num++] =  temp * Nu[ii];
         }
      }
   }

  return;
}



/*
void NurbsSOLID::ShapeFunDerivatives(int* startindex, double* params, double* N, double* dN_dx, double* dN_dy, double* dN_dz, double& Jac)
{
    int ni, nj, nk, ii, jj, kk, loc_num, ind1, ind2, ind3, ROWS = 2;
    double   dR_du[nlbf][3], temp1, temp2, temp3;
    EPOINT *EP;

    ni = startindex[0];
    nj = startindex[1];
    nk = startindex[2];
   
    double** ders1 = new double*[ROWS];
    double** ders2 = new double*[ROWS];
    double** ders3 = new double*[ROWS];

    for(ii=0;ii<ROWS;ii++)
    {
       ders1[ii] = new double[nlbf1];
       ders2[ii] = new double[nlbf2];
       ders3[ii] = new double[nlbf3];
    }

    DersBasisFuns(&(U[0]), U.n, p, params[0], 1, ders1 );
    DersBasisFuns(&(V[0]), V.n, q, params[1], 1, ders2 );
    DersBasisFuns(&(W[0]), W.n, r, params[2], 1, ders3 );

   loc_num = 0;
   for(kk=0;kk<=r;kk++)
   {
      for(jj=0; jj<=q; jj++)
      {
         temp1 = ders2[0][jj] * ders3[0][kk];
         temp2 = ders2[1][jj] * ders3[0][kk];
         temp3 = ders2[0][jj] * ders3[1][kk];

         for(ii=0; ii<=p; ii++)
         {
            //NN[loc_num]        =  ders1[0][ii] * temp1 ;
            dR_du[loc_num][0]  =  ders1[1][ii] * temp1 ;
            dR_du[loc_num][1]  =  ders1[0][ii] * temp2 ;
            dR_du[loc_num][2]  =  ders1[0][ii] * temp3 ;
            loc_num++;
         }
      }
   }

   Matrix3d  dx_du, du_dx;
   dx_du.setZero();
   loc_num = 0;
   for(kk=0;kk<=r;kk++)
   {
      ind1 = nk+kk;
      for(jj=0; jj<=q; jj++)
      {
         ind2 = nj+jj;
         for(ii=0; ii<=p; ii++)
         {
            EP = &(PP[ind1][ni+ii][ind2] );
            
            temp1 = dR_du[loc_num][0];
            temp2 = dR_du[loc_num][1];
            temp3 = dR_du[loc_num][2];

            dx_du(0,0) +=  (EP->x * temp1) ;
            dx_du(1,0) +=  (EP->x * temp2) ;
            dx_du(2,0) +=  (EP->x * temp3) ;
            dx_du(0,1) +=  (EP->y * temp1) ;
            dx_du(1,1) +=  (EP->y * temp2) ;
            dx_du(2,1) +=  (EP->y * temp3) ;
            dx_du(0,2) +=  (EP->z * temp1) ;
            dx_du(1,2) +=  (EP->z * temp2) ;
            dx_du(2,2) +=  (EP->z * temp3) ;

            loc_num++;
         }
      }
   }

   du_dx = dx_du.inverse();

   // Compute derivatives of basis functions w.r.t physical coordinates
   for(loc_num = 0;loc_num<nlbf; loc_num++)
   {
      dN_dx[loc_num] = dR_du[loc_num][0] * du_dx(0,0) + dR_du[loc_num][1] * du_dx(0,1) + dR_du[loc_num][2] * du_dx(0,2);
      dN_dy[loc_num] = dR_du[loc_num][0] * du_dx(1,0) + dR_du[loc_num][1] * du_dx(1,1) + dR_du[loc_num][2] * du_dx(1,2);
      dN_dz[loc_num] = dR_du[loc_num][0] * du_dx(2,0) + dR_du[loc_num][1] * du_dx(2,1) + dR_du[loc_num][2] * du_dx(2,2);
   }

   // only determinant is computed here. The value of the mapping from parametric domain to Master element domain is multiplied in element subroutine

   Jac = dx_du.determinant();

   EP = NULL;

   for(ii=0;ii<ROWS;ii++)
   {
      delete [] ders1[ii];
      delete [] ders2[ii];
      delete [] ders3[ii];
   }
   delete [] ders1;
   delete [] ders2;
   delete [] ders3;

   return;
}
*/



void NurbsSOLID::ShapeFunDerivatives(int* startindex, double* params, double* N, double* dN_dx, double* dN_dy, double* dN_dz, double& Jac)
{

    int ni, nj, nk, ii, jj, kk, loc_num, ind1, ind2, ind3, ROWS = 2;
    double   dR_du[nlbf][3], temp1, temp2, temp3;
    EPOINT *EP;

    ni = startindex[0];
    nj = startindex[1];
    nk = startindex[2];
   
    double** ders1 = new double*[ROWS];
    double** ders2 = new double*[ROWS];
    double** ders3 = new double*[ROWS];

    for(ii=0;ii<ROWS;ii++)
    {
       ders1[ii] = new double[nlbf1];
       ders2[ii] = new double[nlbf2];
       ders3[ii] = new double[nlbf3];
    }

    DersBasisFuns(&(U[0]), U.n, p, params[0], 1, ders1 );
    DersBasisFuns(&(V[0]), V.n, q, params[1], 1, ders2 );
    DersBasisFuns(&(W[0]), W.n, r, params[2], 1, ders3 );

    double  sum1, sum2, sum3, sum_tot, wt, det;

    sum1 = sum2 = sum3 = 0.0;
    sum_tot = 0.0;
    loc_num = 0;

    for(kk=0;kk<=r;kk++)
    {
      ind1 = nk+kk;
      for(jj=0; jj<=q; jj++)
      {
        ind2 = nj+jj;

        temp1 = ders2[0][jj] * ders3[0][kk];
        temp2 = ders2[1][jj] * ders3[0][kk];
        temp3 = ders2[0][jj] * ders3[1][kk];

        for(ii=0; ii<=p; ii++)
        {
          wt = Pw[ind1][ni+ii][ind2].w;

          N[loc_num]        =  ders1[0][ii] * temp1 * wt;

          sum_tot  +=  N[loc_num];
       
          dR_du[loc_num][0]  =  ders1[1][ii] * temp1 * wt;
          dR_du[loc_num][1]  =  ders1[0][ii] * temp2 * wt;
          dR_du[loc_num][2]  =  ders1[0][ii] * temp3 * wt;

          sum1 += dR_du[loc_num][0];
          sum2 += dR_du[loc_num][1];
          sum3 += dR_du[loc_num][2];

          loc_num++;
        }
      }
    }
    
    det = sum_tot*sum_tot;

    for(loc_num = 0;loc_num<nlbf; loc_num++)
    {
      dR_du[loc_num][0] = (dR_du[loc_num][0]*sum_tot - N[loc_num]*sum1)/det;
      dR_du[loc_num][1] = (dR_du[loc_num][1]*sum_tot - N[loc_num]*sum2)/det;
      dR_du[loc_num][2] = (dR_du[loc_num][2]*sum_tot - N[loc_num]*sum3)/det;

      N[loc_num] = N[loc_num]/sum_tot;
    }


    Matrix3d  dx_du, du_dx;
    dx_du.setZero();
    loc_num = 0;
    for(kk=0;kk<=r;kk++)
    {
      ind1 = nk+kk;
      for(jj=0; jj<=q; jj++)
      {
         ind2 = nj+jj;
         for(ii=0; ii<=p; ii++)
         {
            EP = &(PP[ind1][ni+ii][ind2] );
            
            temp1 = dR_du[loc_num][0];
            temp2 = dR_du[loc_num][1];
            temp3 = dR_du[loc_num][2];

            dx_du(0,0) +=  (EP->x * temp1) ;
            dx_du(1,0) +=  (EP->x * temp2) ;
            dx_du(2,0) +=  (EP->x * temp3) ;

            dx_du(0,1) +=  (EP->y * temp1) ;
            dx_du(1,1) +=  (EP->y * temp2) ;
            dx_du(2,1) +=  (EP->y * temp3) ;

            dx_du(0,2) +=  (EP->z * temp1) ;
            dx_du(1,2) +=  (EP->z * temp2) ;
            dx_du(2,2) +=  (EP->z * temp3) ;

            loc_num++;
         }
      }
    }

    du_dx = dx_du.inverse();

    // Compute derivatives of basis functions w.r.t physical coordinates
    for(loc_num = 0;loc_num<nlbf; loc_num++)
    {
      dN_dx[loc_num] = dR_du[loc_num][0] * du_dx(0,0) + dR_du[loc_num][1] * du_dx(0,1) + dR_du[loc_num][2] * du_dx(0,2);
      dN_dy[loc_num] = dR_du[loc_num][0] * du_dx(1,0) + dR_du[loc_num][1] * du_dx(1,1) + dR_du[loc_num][2] * du_dx(1,2);
      dN_dz[loc_num] = dR_du[loc_num][0] * du_dx(2,0) + dR_du[loc_num][1] * du_dx(2,1) + dR_du[loc_num][2] * du_dx(2,2);
    }

    // only determinant is computed here. The value of the mapping from parametric domain to Master element domain is multiplied in element subroutine

    Jac = dx_du.determinant();

    EP = NULL;

   for(ii=0;ii<ROWS;ii++)
   {
      delete [] ders1[ii];
      delete [] ders2[ii];
      delete [] ders3[ii];
   }
   delete [] ders1;
   delete [] ders2;
   delete [] ders3;

   return;
}



//
void NurbsSOLID::deformationGradient(int* startindex, bool flag, double* dN_dx, double* dN_dy, double* dN_dz, double* F, double& detF)
{

    int  ni, nj, nk, ii, jj, kk, loc_num, ind1, ind2;
    ni = startindex[0];
    nj = startindex[1];
    nk = startindex[2];

    EPOINT *EP;
    
    for(ii=0;ii<9;ii++)
      F[ii] = 0.0;
      
    if(flag)
    {
        loc_num = 0;

        for(kk=0; kk<=r; kk++)
        {
           ind1 = nk+kk;
           for(jj=0; jj<=q; jj++)
           {
               ind2 = nj+jj;
               for(ii=0; ii<=p; ii++)
               {
                  EP = &(PP[ind1][ni+ii][ind2] );

                  F[0] += EP->x * dN_dx[loc_num];
                  F[1] += EP->y * dN_dx[loc_num];
                  F[2] += EP->z * dN_dx[loc_num];

                  F[3] += EP->x * dN_dy[loc_num];
                  F[4] += EP->y * dN_dy[loc_num];
                  F[5] += EP->z * dN_dy[loc_num];

                  F[6] += EP->x * dN_dz[loc_num];
                  F[7] += EP->y * dN_dz[loc_num];
                  F[8] += EP->z * dN_dz[loc_num];

                  loc_num++;
               }
           }
        }
        detF = F[0]*(F[4]*F[8] - F[5]*F[7]) - F[3]*(F[1]*F[8] - F[2]*F[7]) + F[6]*(F[1]*F[5] - F[2]*F[4]);
    }
    else //confflag == 2.  calculation of F^-1
    {
        Matrix3d  Finv(3,3), FF(3,3); Finv.setZero();

        loc_num = 0;
        for(kk=0; kk<=r; kk++)
        {
           ind1 = nk+kk;
           for(jj=0; jj<=q; jj++)
           {
               ind2 = nj+jj;
               for(ii=0; ii<=p; ii++)
               {
                  EP = &(PP[ind1][ni+ii][ind2] );

                  Finv(0,0) += EP->x * dN_dx[loc_num];
                  Finv(0,1) += EP->y * dN_dx[loc_num];
                  Finv(0,2) += EP->z * dN_dx[loc_num];

                  Finv(1,0) += EP->x * dN_dy[loc_num];
                  Finv(1,1) += EP->y * dN_dy[loc_num];
                  Finv(1,2) += EP->z * dN_dy[loc_num];

                  Finv(2,0) += EP->x * dN_dz[loc_num];
                  Finv(2,1) += EP->y * dN_dz[loc_num];
                  Finv(2,2) += EP->z * dN_dz[loc_num];

                  loc_num++;
               }
           }
        }
        
        FF = Finv.inverse();

        detF = 1.0 / Finv.determinant();

         F[0] = FF(0,0);         F[3] = FF(1,0);         F[6] = FF(2,0);
         F[1] = FF(0,1);         F[4] = FF(1,1);         F[7] = FF(2,1);
         F[2] = FF(0,2);         F[5] = FF(1,2);         F[8] = FF(2,2);
    }


  return;
}
//


/*
void NurbsSOLID::deformationGradient(int* startindex, bool flag, double* dN_dx, double* dN_dy, double* dN_dz, double* F, double& detF)
{
    int  ni, nj, nk, ii, jj, kk, loc_num, ind1, ind2;
    ni = startindex[0];
    nj = startindex[1];
    nk = startindex[2];
    
    double  dx, dy, dz;

    EPOINT *EP;
    
    for(ii=0;ii<9;ii++)
      F[ii] = 0.0;
    F[0] = F[4] = F[8] = 1.0;

    if(flag)
    {
        loc_num = 0;

        for(kk=0; kk<=r; kk++)
        {
           ind1 = ngbf1m2 * (nk + kk);
           for(jj=0; jj<=q; jj++)
           {
              ind2 = ind1 + ngbf1 * (nj+jj) + ni;
               for(ii=0; ii<=p; ii++)
               {
                  dx = Values[0][ind2+ii];
                  dy = Values[1][ind2+ii];
                  dz = Values[2][ind2+ii];

                  F[0] += dx * dN_dx[loc_num];
                  F[1] += dy * dN_dx[loc_num];
                  F[2] += dz * dN_dx[loc_num];

                  F[3] += dx * dN_dy[loc_num];
                  F[4] += dy * dN_dy[loc_num];
                  F[5] += dz * dN_dy[loc_num];

                  F[6] += dx * dN_dz[loc_num];
                  F[7] += dy * dN_dz[loc_num];
                  F[8] += dz * dN_dz[loc_num];

                  loc_num++;
               }
           }
        }
        detF = F[0]*(F[4]*F[8] - F[5]*F[7]) - F[3]*(F[1]*F[8] - F[2]*F[7]) + F[6]*(F[1]*F[5] - F[2]*F[4]);
    }
    else //confflag == 2.  calculation of F^-1
    {
        Matrix3d  Finv(3,3), FF(3,3); Finv.setZero();

        loc_num = 0;
        for(kk=0; kk<=r; kk++)
        {
           ind1 = nk+kk;
           for(jj=0; jj<=q; jj++)
           {
               ind2 = nj+jj;
               for(ii=0; ii<=p; ii++)
               {
                  EP = &(PP[ind1][ni+ii][ind2] );

                  Finv(0,0) += EP->x * dN_dx[loc_num];
                  Finv(0,1) += EP->y * dN_dx[loc_num];
                  Finv(0,2) += EP->z * dN_dx[loc_num];

                  Finv(1,0) += EP->x * dN_dy[loc_num];
                  Finv(1,1) += EP->y * dN_dy[loc_num];
                  Finv(1,2) += EP->z * dN_dy[loc_num];

                  Finv(2,0) += EP->x * dN_dz[loc_num];
                  Finv(2,1) += EP->y * dN_dz[loc_num];
                  Finv(2,2) += EP->z * dN_dz[loc_num];

                  loc_num++;
               }
           }
        }
        
        FF = Finv.inverse();

        detF = 1.0 / Finv.determinant();

         F[0] = FF(0,0);         F[3] = FF(1,0);         F[6] = FF(2,0);
         F[1] = FF(0,1);         F[4] = FF(1,1);         F[7] = FF(2,1);
         F[2] = FF(0,2);         F[5] = FF(1,2);         F[8] = FF(2,2);
    }


  return;
}
*/


void NurbsSOLID::print2screen()
{
    int index, count, ii, jj, kk;

    for(kk=0;kk<ngbf3;kk++)
    {
        cout << kk << endl;
        for(jj=0;jj<ngbf2;jj++)
        {
            for(ii=0;ii<ngbf1;ii++)
            {
	        PP[kk][ii][jj].print2screen();
            }
            cout << endl;
        }
        cout << endl;
        cout << endl;
    }

   return;
}


void NurbsSOLID::ShapeFunDerivatives3(int face, int* startindex, double* params, double* NN, double* dircos, double& Jac)
{

   assert(face < 7);
   
   int ni, nj, nk, ii, jj, kk, loc_num, ind1, ind2, ind3, ROWS = 2;
   
   double u, v, w, det;
   
   ni = startindex[0];
   nj = startindex[1];
   nk = startindex[2];
   
   u = params[0];
   v = params[1];
   w = params[2];
   
   double** ders1 = new double*[ROWS];
   double** ders2 = new double*[ROWS];

   Matrix2d  dx_du, du_dx;
   EPOINT *EP;

   if(face == 1 || face == 2)
   {
       for(ii=0;ii<ROWS;ii++)
       {
          ders1[ii] = new double[p+1];
          ders2[ii] = new double[q+1];
       }

       DersBasisFuns(&(U[0]), U.n, p, u, 1, ders1 );
       DersBasisFuns(&(V[0]), V.n, q, v, 1, ders2 );

       double   dR_du[(p+1)*(q+1)][2] ;


       loc_num = 0;

       for(jj=0; jj<=q; jj++)
       {
          for(ii=0; ii<=p; ii++)
          {  
             NN[loc_num] = ders1[0][ii] * ders2[0][jj];

             dR_du[loc_num][0] = ders1[1][ii] * ders2[0][jj];
             dR_du[loc_num][1] = ders1[0][ii] * ders2[1][jj];
             loc_num++;
          }
       }

       dx_du.setZero();

       loc_num = 0;
       
       if(face == 1)
         ind3 = 0;
       if(face == 2)
         ind3 = ngbf3-1;

       for(jj=0; jj<=q; jj++)
       {
         ind1 = nj+jj;
         for(ii=0; ii<=p; ii++)
         {
           EP = &(PP[ind3][ni+ii][ind1]);

           dx_du(0,0) +=  (EP->x * dR_du[loc_num][0]) ;
           dx_du(1,0) +=  (EP->x * dR_du[loc_num][1]) ;
           dx_du(0,1) +=  (EP->y * dR_du[loc_num][0]) ;
           dx_du(1,1) +=  (EP->y * dR_du[loc_num][1]) ;

           loc_num++;
         }
       }

       ind1 = ni+p, ind2 = nj+q;

       Jac = dx_du.determinant() * (0.5*(U[ind1+1]-U[ind1])) * (0.5*(V[ind2+1]-V[ind2])) ;
   }
   if(face == 3 || face == 4)
   {
       for(ii=0;ii<ROWS;ii++)
       {
          ders1[ii] = new double[p+1];
          ders2[ii] = new double[r+1];
       }

       DersBasisFuns(&(U[0]), U.n, p, u, 1, ders1 );
       DersBasisFuns(&(W[0]), W.n, r, w, 1, ders2 );

       double   dR_du[(p+1)*(r+1)][2] ;

       loc_num = 0;

       for(kk=0; kk<=r; kk++)
       {
          for(ii=0; ii<=p; ii++)
          {  
             NN[loc_num] = ders1[0][ii] * ders2[0][kk];
             dR_du[loc_num][0] = ders1[1][ii] * ders2[0][kk];
             dR_du[loc_num][1] = ders1[0][ii] * ders2[1][kk];
             loc_num++;
          }
       }

       dx_du.setZero();

       loc_num = 0;
       
       if(face == 3)
         ind3 = 0;
       if(face == 4)
         ind3 = ngbf2-1;

       for(kk=0; kk<=r; kk++)
       {
         ind1 = nk+kk;
         for(ii=0; ii<=p; ii++)
         {
           EP = &(PP[ind1][ni+ii][ind3]);

           dx_du(0,0) +=  (EP->x * dR_du[loc_num][0]) ;
           dx_du(1,0) +=  (EP->x * dR_du[loc_num][1]) ;
           dx_du(0,1) +=  (EP->z * dR_du[loc_num][0]) ;
           dx_du(1,1) +=  (EP->z * dR_du[loc_num][1]) ;

           loc_num++;
         }
       }

       ind1 = ni+p, ind2 = nk+r;

       Jac = dx_du.determinant() * (0.5*(U[ind1+1]-U[ind1])) * (0.5*(W[ind2+1]-W[ind2])) ;
   }
   if(face == 5 || face == 6)
   {
       for(ii=0;ii<ROWS;ii++)
       {
          ders1[ii] = new double[r+1];
          ders2[ii] = new double[q+1];
       }

       DersBasisFuns(&(U[0]), W.n, r, w, 1, ders1 );
       DersBasisFuns(&(V[0]), V.n, q, v, 1, ders2 );

       double   dR_du[(r+1)*(q+1)][2] ;

       loc_num = 0;

       for(jj=0; jj<=q; jj++)
       {
          for(kk=0; kk<=r; kk++)
          {  
             NN[loc_num] = ders1[0][kk] * ders2[0][jj];
             dR_du[loc_num][0] = ders1[1][kk] * ders2[0][jj];
             dR_du[loc_num][1] = ders1[0][kk] * ders2[1][jj];
             loc_num++;
          }
       }

       dx_du.setZero();

       loc_num = 0;
       
       if(face == 5)
         ind3 = 0;
       if(face == 6)
         ind3 = ngbf1-1;

       for(jj=0; jj<=q; jj++)
       {
         ind1 = nj+jj;
         for(kk=0; kk<=r; kk++)
         {
           EP = &(PP[nk+kk][ind3][ind1]);

           dx_du(0,0) +=  (EP->z * dR_du[loc_num][0]) ;
           dx_du(1,0) +=  (EP->z * dR_du[loc_num][1]) ;
           dx_du(0,1) +=  (EP->y * dR_du[loc_num][0]) ;
           dx_du(1,1) +=  (EP->y * dR_du[loc_num][1]) ;

           loc_num++;
         }
       }

       ind1 = nk+r, ind2 = nj+q;

       Jac = dx_du.determinant() * (0.5*(W[ind1+1]-W[ind1])) * (0.5*(V[ind2+1]-V[ind2])) ;
   }   
   
   EP = NULL;

   for(ii=0;ii<ROWS;ii++)
   {
      delete [] ders1[ii];
      delete [] ders2[ii];
   }
   delete [] ders1;
   delete [] ders2;

   return;
}




double  NurbsSOLID::computeValue(int dof, double u, double v, double w)
{
  double  val=0.0, temp;

    double  Nu[nlbf1], Nv[nlbf2], Nw[nlbf3];

    int  ni, nj, nk, row, col, plane, ind1, ind2, kk, ii, jj;

    ni = FindSpan(&(U[0]), U.n, p, u) - p;
    nj = FindSpan(&(V[0]), V.n, q, v) - q;
    nk = FindSpan(&(W[0]), W.n, r, w) - r;

    BasisFuns(&(U[0]), U.n, p, u, Nu) ;
    BasisFuns(&(V[0]), V.n, q, v, Nv) ;
    BasisFuns(&(W[0]), W.n, r, w, Nw) ;

    val = 0.0;
    dof -= 1;

    for(kk=0;kk<=r;kk++)
    {
      ind1 = ngbf1m2 * (nk + kk);
      for(jj=0;jj<=q;jj++)
      {
        ind2 = ind1 + ngbf1 * (nj+jj) + ni;
        temp = Nw[kk] * Nv[jj];
        for(ii=0;ii<=p;ii++)
        {
          val += Nu[ii] * temp * Values[dof][ind2+ii];
        }
      }
    }

  return val;
}





double  NurbsSOLID::computeValueAndShanpeFns(int dof, double u, double v, double w, double* NN)
{
  double  val=0.0, temp;
  
    double  Nu[nlbf1], Nv[nlbf2], Nw[nlbf3], temp1, temp2;

    int  ni, nj, nk, row, col, plane, ind1, ind2, kk, ii, jj, count;

    ni = FindSpan(&(U[0]), U.n, p, u) - p;
    nj = FindSpan(&(V[0]), V.n, q, v) - q;
    nk = FindSpan(&(W[0]), W.n, r, w) - r;
     
    BasisFuns(&(U[0]), U.n, p, u, Nu) ;
    BasisFuns(&(V[0]), V.n, q, v, Nv) ;
    BasisFuns(&(W[0]), W.n, r, w, Nw) ;

    val = 0.0;
    count=0;

    dof -= 1;
    for(kk=0;kk<=r;kk++)
    {
      ind1 = ngbf1m2 * (nk + kk);
      for(jj=0;jj<=q;jj++)
      {
        ind2 = ind1 + ngbf1 * (nj+jj) + ni;
        temp1 = Nw[kk] * Nv[jj];
        for(ii=0;ii<=p;ii++)
        {
          temp2 = temp1 *  Nu[ii];
          NN[count++] =  temp2;
          val += temp2 * Values[dof][ind2+ii];
        }
      }
    }

  return val;
}





