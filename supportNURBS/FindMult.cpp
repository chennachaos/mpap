#include <iostream>
#include <math.h>
#include "NURBS_1D.h"

using namespace std;



//  Determines the knot multiplicity for a given value of parameter 'u'

//int FindMult(KNOTVECTOR& U, DEGREE deg, double ubar)


int FindMult(double* U, int Un, DEGREE deg, double ubar)
{
  int m = Un-1;
  int n = m-deg-1;

  if(CompareDoubles(ubar, U[0]) || CompareDoubles(ubar, U[m]))  // special case
    return deg+1 ;

  int s = 0;
  int span = FindSpan(U, Un, deg, ubar);
  
  for(int i=0;i<deg;i++)
  {
    if(CompareDoubles(ubar, U[span-i]))  
      s++ ;
    else  
      s = s;
  }
  return s ;
}
