#include <iostream>
#include <math.h>
#include "NURBS_1D.h"

using namespace std;


//   Determines if a value "u" is a knot of KNOTVECTOR U

bool IsKnot(double* U, int Un, double u)
{
  int ii, count = 0;

  for(ii=0;ii<Un;ii++)
  {
    if(CompareDoubles(u, U[ii]))
      return true ;
  }

  return false;
}

