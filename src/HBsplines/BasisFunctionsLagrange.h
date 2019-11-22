
#ifndef incl_BasisFunctionsLagrange_h
#define incl_BasisFunctionsLagrange_h


/**
  Computes function values, first derivative and second derivative of
  univariate Lagrange polynomials of a given order p for
*/
void Lagrange_BasisFuns1D(int p, double xi, double* N, double* dN_dxi, double* d2N_dxi2);

/**
  Computes function values and first derivative of
  univariate Lagrange polynomials of a given order p
*/
void Lagrange_BasisFuns1D(int p, double xi, double* N, double* dN_dxi);

/**
  Computes function values of
  univariate Lagrange polynomials of a given order p
*/
void Lagrange_BasisFuns1D(int p, double xi, double* N);


void computeLagrangeBFsLine1D(int p, double uu, double *xx, double *N, double *dN_dx, double& Jac);

void computeLagrangeBFsLine2D(int p, double uu, double *xx, double* yy, double *N, double *dN_dx, double& Jac);

void computeLagrangeBFsLine3D(int p, double uu, double *xx, double* yy, double* zz, double *N, double *dN_dx, double& Jac);

/**
  Computes function values and first derivative of
  Lagrange polynomials of a given order p for Triangular elements
*/
void LagrangeBasisFunsTria(int p, double xi, double zeta, double* N, double* dN_dxi, double* dN_dzeta);

/**
  Computes function values and first derivative of
  Lagrange polynomials of a given order p for Quadrilateral elements
*/
void LagrangeBasisFunsQuad(int p, double xi, double zeta, double* N, double* dN_dxi, double* dN_dzeta);

/**
  Computes function values and first derivative of
  Lagrange polynomials of a given order p for Tetrahedral elements
*/
void LagrangeBasisFunsTet(int p, double xi, double zeta, double eta, double* N, double* dN_dxi, double* dN_dzeta, double* dN_deta);

/**
  Computes function values and first derivative of
  Lagrange polynomials of a given order p for Hexahedral elements
*/
void LagrangeBasisFunsHex(int p, double xi, double zeta, double eta, double* N, double* dN_dxi, double* dN_dzeta, double* dN_deta);





#endif


