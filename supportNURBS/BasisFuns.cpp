#include <iostream>
#include <math.h>
#include "NURBS_1D.h"

using namespace std;


/*
   Computes the non-zero basis functions into N of size curve_deg+1
   Algorithm A2.2 on pg#70 of the Nurbs Book
   u  the parametric value
   i  the non-zero span of the basis functions
   N  the non-zero basis functions
*/


//void BasisFuns(KNOTVECTOR& U, DEGREE p, double u, VectorArray<double>& N)


void BasisFuns(double* U, int Un, DEGREE p, double u, double* N)
{
  vector<double>  left(p+1), right(p+1);
  double temp, saved ;

  int span, j, r;
  
  span = FindSpan(U, Un, p, u);

  N[0] = 1.0 ;
  
  for(j=1;j<=p;j++)
  {
     left[j] = u-U[span+1-j] ;
     right[j] = U[span+j]-u ;
     saved = 0.0 ;
     for(r=0;r<j;r++)
     {
        temp = N[r]/(right[r+1]+left[j-r]) ;
        N[r] = saved+right[r+1] * temp ;
        saved = left[j-r] * temp ;
     }
     N[j] = saved ;
  }
  
  return;
}
