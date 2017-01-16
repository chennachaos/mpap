#include <iostream>
#include <math.h>
#include "NURBS_1D.h"

using namespace std;



/*
  Computes the k-th derivatives of Basis function "i" at u of the NURBS curve
  Using the recursive definition of derivatives equation 2.9 on pg#61 of the NURBS book.
  i     the number of NURBS basis function for which to calculate the derivatives
  u     the parametric value
  k     the degree of the derivation
  ders  A vector containing the derivatives(0 through k) of the curve.
*/


double DersOneBasisFun_recurs(KNOTVECTOR& U, DEGREE deg, int i, double u, int der_order)
{
  double der_value = 0.0;
  if (der_order==0)
  {
    der_value = OneBasisFun_recurs(U, deg, i, u);
    return der_value;
  }
  else
  {
    double dummy1 = 0.0, dummy2 = 0.0, dummy3 = 0.0, dummy4 = 0.0;

    dummy1 = DersOneBasisFun_recurs(U, deg-1, i, u, der_order-1);
    dummy2 = DersOneBasisFun_recurs(U, deg-1, i+1, u, der_order-1);

    if (CompareDoubles(dummy1, 0.0))
      dummy3 = 0.0;
    else
      dummy3 = dummy1/(U[i+deg]-U[i]);

    if (CompareDoubles(dummy2, 0.0))
      dummy4 = 0.0;
    else
      dummy4 = dummy2/(U[i+deg+1]-U[i+1]);

    der_value = deg*(dummy3-dummy4);
    return der_value;
  }
}
