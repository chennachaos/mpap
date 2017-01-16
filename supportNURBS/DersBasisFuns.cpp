#include <iostream>
#include <math.h>
#include "NURBS_1D.h"

using namespace std;



/*
  Computes the derivatives of Non-zero Basis functions at u of the NURBS curve
  Algorithm A2.3 on pg#72 of the NURBS book.

  n     the degree of the derivation
  u     the parametric value
  span  the span for the basis functions
  ders  A matrix containing the derivatives of the curve.
*/


void DersBasisFuns(double* U, int Un, DEGREE p, double u, int n, double** ders)
{
  int j, r, s1, s2, rk, pk, j1, j2, k, span;

  //double left[p+1], right[p+1];
  //double ndu[p+1][p+1], a[p+1][p+1] ;
  
  vector<double>  left(p + 1), right(p + 1);
  vector<vector<double> >  ndu, a;

  k = p + 1;

  ndu.resize(k);
  a.resize(k);

  for(j=0; j<k; j++)
  {
    ndu[j].resize(k);
    a[j].resize(k);
  }
  
  double saved, temp, d;
  
  span = FindSpan(U, Un, p, u);

  ndu[0][0]=1.0 ;
  for(j=1; j<= p; j++)
  {
    left[j] = u-U[span+1-j] ;
    right[j] = U[span+j]-u ;
    saved = 0.0 ;
    
    for(r=0;r<j;r++)
    {
      // Lower triangle
      ndu[j][r] = right[r+1] + left[j-r] ;
      temp = ndu[r][j-1]/ndu[j][r] ;

      // Upper triangle
      ndu[r][j] = saved + right[r+1] * temp ;
      saved = left[j-r] * temp ;
    }

    ndu[j][j] = saved ;
  }

  for(j=0;j<=p;j++)
    ders[0][j] = ndu[j][p];

  // Compute the derivatives
 
  
  for(r=0;r<=p;r++)
  {
    s1=0; s2=1 ; // alternate rows in array a
    a[0][0] = 1.0 ;
    
    // Compute the kth derivative
    for(k=1;k<=n;k++)
    {
      d=0.0 ;
      rk = r-k ;
      pk = p-k ;

      if(r>=k)
      {
	  a[s2][0] = a[s1][0]/ndu[pk+1][rk] ;
	  d = a[s2][0] * ndu[rk][pk] ;
      }

      if(rk >= -1)
        j1 = 1;
      else
        j1 = -rk;

      if(r-1 <= pk)
        j2 = k-1 ;
      else
        j2 = p-r ;

      for(j=j1;j<=j2;j++)
      {
	 a[s2][j] = (a[s1][j]-a[s1][j-1])/ndu[pk+1][rk+j];
	 d += a[s2][j] * ndu[rk+j][pk];
      }
      
      if(r<=pk)
      {
        a[s2][k] = -a[s1][k-1]/ndu[pk+1][r] ;
        d += a[s2][k] * ndu[r][pk] ;
      }
      ders[k][r] = d ;
      j = s1 ; s1 = s2 ; s2 = j ; // Switch rows
    }
  }
  // Multiply through by the correct factors
  r = p ;
  for(k=1;k<=n;k++)
  {
    for(j=0;j<=p;j++)
      ders[k][j] *= r ;
    r *= (p - k) ;
  }

 return;

}
