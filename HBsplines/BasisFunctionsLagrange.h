
#ifndef incl_BasisFunctionsLagrange_h
#define incl_BasisFunctionsLagrange_h


//#include "headersBasic.h"



void Lagrange_BasisFuns1D(int p, double xi, double* N, double* dN_dxi, double* d2N_dxi2);

void Lagrange_BasisFuns1D(int p, double xi, double* N, double* dN_dxi);

void Lagrange_BasisFuns1D(int p, double xi, double* N);

void computeLagrangeBFsLine1D(int p, double uu, double *xx, double *N, double *dN_dx, double& Jac);

void computeLagrangeBFsLine2D(int p, double uu, double *xx, double* yy, double *N, double *dN_dx, double& Jac);

void computeLagrangeBFsLine3D(int p, double uu, double *xx, double* yy, double* zz, double *N, double *dN_dx, double& Jac);

void  computeLagrangeBFs3(double uu, double *xx, double *yy, double *N);

void LagrangeBasisFunsTria(int p, double xi, double zeta, double* N, double* dN_dxi, double* dN_dzeta);

void LagrangeBasisFunsQuad(int p, double xi, double zeta, double* N, double* dN_dxi, double* dN_dzeta);

void LagrangeBasisFunsTet(int p, double xi, double zeta, double eta, double* N, double* dN_dxi, double* dN_dzeta, double* dN_deta);

void LagrangeBasisFunsHex(int p, double xi, double zeta, double eta, double* N, double* dN_dxi, double* dN_dzeta, double* dN_deta);





#endif


