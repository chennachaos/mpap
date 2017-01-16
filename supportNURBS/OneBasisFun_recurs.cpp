#include <iostream>
#include <math.h>
#include "NURBS_1D.h"

using namespace std;



/*
  Computes the i'th Basis function for the knot vector U, degree p of the curve at parameter u 
  The basis function is noted as Nip(u)
  Recursive function based on equation 2.5 in pg#50 of the Nurbs Book
  u  the parametric variable
  i  specifies which basis function to compute
  p  the degree to which the basis function is computed
  return the value of  Nip
*/


double OneBasisFun_recurs(KNOTVECTOR& U, DEGREE deg, int i, double u)
{
  int m = U.n-1;
  int n = m-deg-1;
  double Nip = 0.0;

  if (i == n && CompareDoubles(u, U[m]))  // Special Cases
  {
    Nip = 1.0;
    return Nip;
  }

  if (CompareDoubles(deg, 0.0))
  {
     if ( u >= U[i] && u < U[i+1] )
     {
        Nip = 1.0;
	return Nip;
     }
     else
     {
       Nip = 0.0 ; 
       return Nip;
     }
  }
  else
  {
    double dummy1 = 0.0, dummy2 = 0.0, dummy3 = 0.0, dummy4 = 0.0;

    dummy1 = OneBasisFun_recurs(U, deg-1, i, u);
    dummy2 = OneBasisFun_recurs(U, deg-1, i+1, u);

    if (CompareDoubles(dummy1, 0.0))
      dummy3=0.0;
    else
      dummy3=dummy1*(u-U[i])/(U[i+deg]-U[i]);

    if (CompareDoubles(dummy2, 0.0))
      dummy4=0.0;
    else
      dummy4=dummy2*(U[i+deg+1]-u)/(U[i+deg+1]-U[i+1]);

    Nip = dummy3 + dummy4;
    return Nip;
  }
}
