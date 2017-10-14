
#include "BasisFunctionsLagrange.h"

#include <math.h>
#include <cmath>
#include <stdio.h>
#include <vector>

using std::vector;


void Lagrange_BasisFuns1D(int p, double xi, double* N, double* dN_dxi, double* d2N_dxi2)
{
  double  fact1, fact2, val1, val2, val3, val4;

  switch(p)
  {
      case 0:

          N[0] = 1.0;

          dN_dxi[0] = 0.0;
    
          d2N_dxi2[0] = 0.0;

      break;

      case 1:

          N[0] = 0.5*(1.0 - xi);
          N[1] = 0.5*(1.0 + xi);
    
          dN_dxi[0] = -0.5;
          dN_dxi[1] =  0.5;
    
          d2N_dxi2[0] = 0.0;
          d2N_dxi2[1] = 0.0;

      break;

      case 2:

          val1 = xi*xi;

          N[0] = -0.5 * (xi - val1);
          N[1] = 1.0 - val1;
          N[2] = 0.5 *(xi + val1);

          val1 = 2.0*xi;

          dN_dxi[0] = -0.5*(1.0 - val1);
          dN_dxi[1] = -val1;
          dN_dxi[2] =  0.5*(1.0 + val1);

          d2N_dxi2[0] =  1.0;
          d2N_dxi2[1] = -2.0;
          d2N_dxi2[2] =  1.0;

      break;

      case 3:

          fact1 = 9.0/16.0;
          fact2 = 27.0/16.0;
          val1  = xi*xi;

          N[0] = -fact1 * (1 - xi)   * (1/9 - val1);
          N[1] =  fact2 * (1 - val1) * (1/3 - xi);
          N[2] =  fact2 * (1 - val1) * (1/3 + xi);
          N[3] = -fact1 * (1 + xi)   * ( 1/9 - val1);

          val2 = 3.0*val1;

          dN_dxi[0] = -fact1*(-1/9 - 2.0*xi   +  val2);
          dN_dxi[1] =  fact2*(-1   - 2.0/3*xi +  val2);
          dN_dxi[2] =  fact2*(1    - 2.0/3*xi -  val2);
          dN_dxi[3] = -fact1*(1/9  - 2.0*xi   -  val2);
    
          val2 = 6.0*xi;
    
          d2N_dxi2[0] = -fact1 * (-2   + val2);
          d2N_dxi2[1] =  fact2 * (-2/3 + val2);
          d2N_dxi2[2] =  fact2 * (-2/3 - val2);
          d2N_dxi2[3] = -fact1 * (-2   - val2);

      break;

      case 4:
     
          fact1 = 2.0/3.0;
          fact2 = 8.0/3.0;
          val1 = xi*xi;
          val2 = val1*xi;
          val3 = val2*xi;
    
          N[0] =  fact1 * (0.25*xi  - 0.25*val1  -  val2     + val3);
          N[1] = -fact2 * (0.5 *xi  - val1       -  0.5*val2 + val3);
          N[2] =    4.0 * (0.25     - 1.25*val1  -  0        + val3);
          N[3] =  fact2 * (0.5*xi   + val1       -  0.5*val2 - val3);
          N[4] = -fact1 * (0.25*xi  + 0.25*val1  -  val2     - val3);
    
          val4 = 4.0*val2;

          dN_dxi[0] =  fact1 * (0.25 - 0.5*xi  - 3.0*val1  + val4);
          dN_dxi[1] = -fact2 * (0.5  - 2.0*xi  - 1.5*val1  + val4);
          dN_dxi[2] =    4.0 * (  0  - 2.5*xi  -   0.0     + val4);
          dN_dxi[3] =  fact2 * (0.5  + 2.0*xi  - 1.5*val1  - val4);
          dN_dxi[4] = -fact1 * (0.25 + 0.5*xi  - 3.0*val1  - val4);
    
          val4 = 12.0*val1;
    
          d2N_dxi2[0] =  fact1 * (-0.5  -  6.0*xi  +  val4);
          d2N_dxi2[1] = -fact2 * (-2.0  -  3.0*xi  +  val4);
          d2N_dxi2[2] =    4.0 * (-2.5  -  0.0     +  val4);
          d2N_dxi2[3] =  fact2 * ( 2.0  -  3.0*xi  -  val4);
          d2N_dxi2[4] = -fact1 * ( 0.5  -  6.0*xi  -  val4);

      break;

      default:

          printf("no basis functions defined for this degree = %5d \n", p);

      break;
  }

   return;
}



void Lagrange_BasisFuns1D(int p, double xi, double* N, double* dN_dxi)
{
  double  fact1, fact2, val1, val2, val3, val4;

  switch(p)
  {
      case 0:

          N[0] = 1.0;

          dN_dxi[0] = 0.0;

      break;

      case 1:

          N[0] = 0.5*(1.0 - xi);
          N[1] = 0.5*(1.0 + xi);
    
          dN_dxi[0] = -0.5;
          dN_dxi[1] =  0.5;

      break;

      case 2:

          val1 = xi*xi;

          N[0] = -0.5 * (xi - val1);
          N[1] = 1.0 - val1;
          N[2] = 0.5 *(xi + val1);
    
          val1 = 2.0*xi;
    
          dN_dxi[0] = -0.5*(1.0 - val1);
          dN_dxi[1] = -val1;
          dN_dxi[2] =  0.5*(1.0 + val1);

      break;

      case 3:

          fact1 = 9.0/16.0;
          fact2 = 27.0/16.0;
          val1 = xi*xi;
    
          N[0] = -fact1 * (1 - xi)   * (1/9 - val1);
          N[1] =  fact2 * (1 - val1) * (1/3 - xi);
          N[2] =  fact2 * (1 - val1) * (1/3 + xi);
          N[3] = -fact1 * (1 + xi)   * ( 1/9 - val1);
    
          val2 = 3.0*val1;

          dN_dxi[0] = -fact1*(-1/9 - 2.0*xi   +  val2);
          dN_dxi[1] =  fact2*(-1   - 2.0/3*xi +  val2);
          dN_dxi[2] =  fact2*(1    - 2.0/3*xi -  val2);
          dN_dxi[3] = -fact1*(1/9  - 2.0*xi   -  val2);

      break;

      case 4:

          fact1 = 2.0/3.0;
          fact2 = 8.0/3.0;
          val1 = xi*xi;
          val2 = val1*xi;
          val3 = val2*xi;
    
          N[0] =  fact1 * (0.25*xi  - 0.25*val1  -  val2     + val3);
          N[1] = -fact2 * (0.5 *xi  - val1       -  0.5*val2 + val3);
          N[2] =    4.0 * (0.25     - 1.25*val1  -  0        + val3);
          N[3] =  fact2 * (0.5*xi   + val1       -  0.5*val2 - val3);
          N[4] = -fact1 * (0.25*xi  + 0.25*val1  -  val2     - val3);
    
          val4 = 4.0*val2;

          dN_dxi[0] =  fact1 * (0.25 - 0.5*xi  - 3.0*val1  + val4);
          dN_dxi[1] = -fact2 * (0.5  - 2.0*xi  - 1.5*val1  + val4);
          dN_dxi[2] =    4.0 * (  0  - 2.5*xi  -   0.0     + val4);
          dN_dxi[3] =  fact2 * (0.5  + 2.0*xi  - 1.5*val1  - val4);
          dN_dxi[4] = -fact1 * (0.25 + 0.5*xi  - 3.0*val1  - val4);

      break;

      default:

          printf("no basis functions defined for this degree = %5d \n", p);

      break;
  }

  return;
}



void Lagrange_BasisFuns1D(int p, double xi, double* N)
{
  double  fact1, fact2, val1, val2, val3, val4;

  switch(p)
  {
      case 0:

          N[0] = 1.0;

      break;

      case 1:

          N[0] = 0.5*(1.0 - xi);
          N[1] = 0.5*(1.0 + xi);

      break;

      case 2:

          val1 = xi*xi;

          N[0] = -0.5 * (xi - val1);
          N[1] = 1.0 - val1;
          N[2] = 0.5 *(xi + val1);

      break;

      case 3:

          fact1 = 9.0/16.0;
          fact2 = 27.0/16.0;
          val1 = xi*xi;
    
          N[0] = -fact1 * (1 - xi)   * (1/9 - val1);
          N[1] =  fact2 * (1 - val1) * (1/3 - xi);
          N[2] =  fact2 * (1 - val1) * (1/3 + xi);
          N[3] = -fact1 * (1 + xi)   * ( 1/9 - val1);
      break;

      case 4:

          fact1 = 2.0/3.0;
          fact2 = 8.0/3.0;
          val1 = xi*xi;
          val2 = val1*xi;
          val3 = val2*xi;
    
          N[0] =  fact1 * (0.25*xi  - 0.25*val1  -  val2     + val3);
          N[1] = -fact2 * (0.5 *xi  - val1       -  0.5*val2 + val3);
          N[2] =    4.0 * (0.25     - 1.25*val1  -  0        + val3);
          N[3] =  fact2 * (0.5*xi   + val1       -  0.5*val2 - val3);
          N[4] = -fact1 * (0.25*xi  + 0.25*val1  -  val2     - val3);

      break;

      default:

          printf("no basis functions defined for this degree = %5d \n", p);

      break;
  }

   return;
}



void computeLagrangeBFsLine1D(int p, double uu, double *xx, double *N, double *dN_dx, double& Jac)
{
  int  ii, jj, count, nlbf = p+1;

  double dx_du ;
  vector<double>  dN1(nlbf);

  Lagrange_BasisFuns1D(p, uu, N, &dN1[0]);

  if(p == 0)
  {
    Jac = 1.0;
  }
  else
  {
    Jac = 0.0;
    for(ii=0; ii<nlbf; ii++)
      Jac +=  (xx[ii] * dN1[ii]);
  }
  
  dx_du = 1.0/Jac ;

  // Compute derivatives of basis functions w.r.t physical coordinates
  for(ii=0; ii<nlbf; ii++)
    dN_dx[ii] = dN1[ii] * dx_du ;
  
  return;
}



void computeLagrangeBFsLine2D(int p, double uu, double *xx, double* yy, double *N, double *dN_dx, double& Jac)
{
   int  ii, jj, count, nlbf = p+1;

   double  du_dx, dx, dy ;
   vector<double>  dN1(nlbf);

   Lagrange_BasisFuns1D(p, uu, N, &dN1[0]);

   if(p == 0)
   {
     Jac = 1.0;
   }
   else
   {
     dx = dy = 0.0;
     for(ii=0; ii<nlbf; ii++)
     {
       dx +=  (xx[ii] * dN1[ii]);
       dy +=  (yy[ii] * dN1[ii]);
     }
     Jac = sqrt(dx*dx+dy*dy);
   }

   du_dx = 1.0/Jac ;

  // Compute derivatives of basis functions w.r.t physical coordinates
  for(ii=0; ii<nlbf; ii++)
    dN_dx[ii] = dN1[ii] * du_dx ;
  
  return;
}



void computeLagrangeBFsLine3D(int p, double uu, double *xx, double* yy, double* zz, double *N, double *dN_dx, double& Jac)
{
   int  ii, jj, count, nlbf = p+1;

   double  du_dx, dx, dy, dz ;
   vector<double>  dN1(nlbf);

   Lagrange_BasisFuns1D(p, uu, N, &dN1[0]);

   if(p == 0)
   {
     Jac = 1.0;
   }
   else
   {
     dx = dy = dz = 0.0;
     for(ii=0; ii<nlbf; ii++)
     {
       dx +=  (xx[ii] * dN1[ii]);
       dy +=  (yy[ii] * dN1[ii]);
       dz +=  (zz[ii] * dN1[ii]);
     }
     Jac = sqrt(dx*dx+dy*dy+dz*dz);
   }

   du_dx = 1.0/Jac ;

  // Compute derivatives of basis functions w.r.t physical coordinates
  for(ii=0; ii<nlbf; ii++)
    dN_dx[ii] = dN1[ii] * du_dx ;
  
  return;
}




void LagrangeBasisFunsTria(int p, double xi, double zeta, double* N, double* dN_dxi, double* dN_dzeta)
{
  double  fact1, fact2, val1, val2, val3, val4;
  
  switch(p)
  {
      case 0:

          N[0] = 1.0;

          dN_dxi[0] = 0.0;

      break;

      case 1:

          N[0] = 1.0 - xi - zeta;
          N[1] = xi;
          N[2] = zeta;

          dN_dxi[0] = -1.0;
          dN_dxi[1] =  1.0;
          dN_dxi[2] =  0.0;

          dN_dzeta[0] = -1.0;
          dN_dzeta[1] =  0.0;
          dN_dzeta[2] =  1.0;

      break;

      case 2:

          val1 = xi*xi;

          N[0] = -0.5 * (xi - val1);
          N[1] = 1.0 - val1;
          N[2] = 0.5 *(xi + val1);
          N[3] = -0.5 * (xi - val1);
          N[4] = 1.0 - val1;
          N[5] = 0.5 *(xi + val1);

    
          val1 = 2.0*xi;
    
          dN_dxi[0] = -0.5*(1.0 - val1);
          dN_dxi[1] = -val1;
          dN_dxi[2] =  0.5*(1.0 + val1);
          dN_dxi[3] = -1.0;
          dN_dxi[4] =  1.0;
          dN_dxi[5] =  0.0;

          dN_dzeta[0] = -1.0;
          dN_dzeta[1] =  0.0;
          dN_dzeta[2] =  1.0;
          dN_dzeta[3] = -1.0;
          dN_dzeta[4] =  0.0;
          dN_dzeta[5] =  1.0;

      break;

      default:

          printf("no basis functions defined for this degree = %5d \n", p);

      break;
  }

}



void LagrangeBasisFunsQuad(int p, double xi, double eta, double* N, double* dN_dxi, double* dN_deta)
{
  double  fact1, fact2, v1, v2, v3, v4, v5, v6;

  switch(p)
  {
      case 0:

          N[0] = 1.0;

          dN_dxi[0] = 0.0;
          dN_deta[0] = 0.0;

      break;

      case 1:

          v1 = 1.0 - xi;
          v2 = 1.0 + xi;
          v3 = 1.0 - eta;
          v4 = 1.0 + eta;

          N[0] = 0.25*v1*v3;
          N[1] = 0.25*v2*v3;
          N[2] = 0.25*v1*v4;
          N[3] = 0.25*v2*v4;

          dN_dxi[0] = -0.25*v3;
          dN_dxi[1] =  0.25*v3;
          dN_dxi[2] = -0.25*v4;
          dN_dxi[3] =  0.25*v4;

          dN_deta[0] = -0.25*v1;
          dN_deta[1] = -0.25*v2;
          dN_deta[2] =  0.25*v1;
          dN_deta[3] =  0.25*v2;

      break;

      case 2:
          
          fact1 = xi*xi;
          fact2 = eta*eta;

          v1 = 0.5*(fact1 - xi);
          v2 = 1.0 - fact1;
          v3 = 0.5*(fact1 + xi);

          v4 = 0.5*(fact2 - eta);
          v5 = 1.0 - fact2;
          v6 = 0.5*(fact2 + eta);

          N[0] = v4*v1;
          N[1] = v4*v2;
          N[2] = v4*v3;
          N[3] = v5*v1;
          N[4] = v5*v2;
          N[5] = v5*v3;
          N[6] = v6*v1;
          N[7] = v6*v2;
          N[8] = v6*v3;

          dN_dxi[0] = -0.25*v3;
          dN_dxi[1] =  0.25*v3;
          dN_dxi[2] = -0.25*v4;
          dN_dxi[3] =  0.25*v4;

          dN_deta[0] = -0.25*v1;
          dN_deta[1] = -0.25*v2;
          dN_deta[2] =  0.25*v1;
          dN_deta[3] =  0.25*v2;

      break;

      default:

          printf("no basis functions defined for this degree = %5d \n", p);

      break;
  }

  return;
}




void LagrangeBasisFunsTet(int p, double xi1, double xi2, double xi3, double* N, double* dN_dxi1, double* dN_dxi2, double* dN_dxi3)
{
  double  fact1, fact2, val1, val2, val3, val4;
  
  switch(p)
  {
      case 0:

          N[0] = 1.0;

          dN_dxi1[0] = 0.0;
          dN_dxi2[0] = 0.0;
          dN_dxi3[0] = 0.0;

      break;

      case 1:

          N[0] = xi1;
          N[1] = xi3;
          N[2] = xi2;
          N[3] = 1.0 - xi1 - xi2 - xi3;

          dN_dxi1[0] =  1.0;
          dN_dxi1[1] =  0.0;
          dN_dxi1[2] =  0.0;
          dN_dxi1[3] = -1.0;

          dN_dxi2[0] =  0.0;
          dN_dxi2[1] =  0.0;
          dN_dxi2[2] =  1.0;
          dN_dxi2[3] = -1.0;

          dN_dxi3[0] =  0.0;
          dN_dxi3[1] =  1.0;
          dN_dxi3[2] =  0.0;
          dN_dxi3[3] = -1.0;

      break;

      default:

          printf("no basis functions defined for this degree = %5d \n", p);

      break;
  }

}



void LagrangeBasisFunsHex(int p, double xi1, double xi2, double xi3, double* N, double* dN_dxi1, double* dN_dxi2, double* dN_dxi3)
{
  double  v11, v12, v21, v22, v31, v32;

  switch(p)
  {
      case 0:

          N[0] = 1.0;

          dN_dxi1[0] = 0.0;
          dN_dxi2[0] = 0.0;
          dN_dxi3[0] = 0.0;

      break;

      case 1:

          v11 = 1.0 - xi1;
          v12 = 1.0 + xi1;
          v21 = 1.0 - xi2;
          v22 = 1.0 + xi2;
          v31 = 1.0 - xi3;
          v32 = 1.0 + xi3;

          N[0] = 0.125*v11*v21*v31;
          N[1] = 0.125*v12*v21*v31;
          N[2] = 0.125*v11*v22*v31;
          N[3] = 0.125*v12*v22*v31;
          N[4] = 0.125*v11*v21*v32;
          N[5] = 0.125*v12*v21*v32;
          N[6] = 0.125*v11*v22*v32;
          N[7] = 0.125*v12*v22*v32;

          dN_dxi1[0] = -0.125*v21*v31;
          dN_dxi1[1] =  0.125*v21*v31;
          dN_dxi1[2] = -0.125*v22*v31;
          dN_dxi1[3] =  0.125*v22*v31;
          dN_dxi1[4] = -0.125*v21*v32;
          dN_dxi1[5] =  0.125*v21*v32;
          dN_dxi1[6] = -0.125*v22*v32;
          dN_dxi1[7] =  0.125*v22*v32;

          dN_dxi2[0] = -0.125*v11*v31;
          dN_dxi2[1] = -0.125*v12*v31;
          dN_dxi2[2] =  0.125*v11*v31;
          dN_dxi2[3] =  0.125*v12*v31;
          dN_dxi2[4] = -0.125*v11*v32;
          dN_dxi2[5] = -0.125*v12*v32;
          dN_dxi2[6] =  0.125*v11*v32;
          dN_dxi2[7] =  0.125*v12*v32;

          dN_dxi3[0] = -0.125*v11*v21;
          dN_dxi3[1] = -0.125*v12*v21;
          dN_dxi3[2] = -0.125*v11*v22;
          dN_dxi3[3] = -0.125*v12*v22;
          dN_dxi3[4] =  0.125*v11*v21;
          dN_dxi3[5] =  0.125*v12*v21;
          dN_dxi3[6] =  0.125*v11*v22;
          dN_dxi3[7] =  0.125*v12*v22;

      break;

      default:

          printf("no basis functions defined for this degree = %5d \n", p);

      break;
  }

  return;
}









