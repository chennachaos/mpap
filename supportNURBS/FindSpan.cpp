#include <iostream>
#include <math.h>
#include "NURBS_1D.h"

using namespace std;

/*
  Determines the knot span index for a given value of parameter 'u'
  Algorithm A2.1 on pg#68 of the Nurbs Book
*/

// Un = size of array U

int FindSpan(double* U, int Un, int deg, double u)
{
//  cout << Un << '\t' << deg << '\t' << u << endl;

  int n = Un-deg-2;

  if(CompareDoubles(u, U[n+1]))  // special case
    return n ;

  int low = deg ;
  int high = n+1 ;
  int mid = 0.5*(low+high) ;

  while(u<U[mid] || u>=U[mid+1])
  {
    if(u<U[mid])  high = mid ;
    else  low = mid ;
    mid = 0.5*(low+high) ;
  }

  return mid ;
}
