
#include "BasisFunctionsBSpline.h"

#include <vector>

using std::vector;

//////////////////////////////////////////////////////////////
// compute univariate B-Spline basis functions for 
// at parameter 'u'
// degree 'p'
// starting at knot 'start'
// with a knot span of 'incr'
// and store in array 'N'
//////////////////////////////////////////////////////////////


void HB_BasisFuns(int p, double start, double incr, double u, double* N)
{
  int  r, j(p + 1);

  double temp, saved;
  
  std::vector<double>  left(j), right(j);
  
  N[0] = 1.0 ;
  for(j=1;j<=p;j++)
  {
     left[j] = u - (start-incr*(j-1)) ;
     right[j] = start + incr*j-u ;
     saved = 0.0 ;
     for(r=0;r<j;r++)
     {
        temp = N[r]/(right[r+1]+left[j-r]) ;
        N[r] = saved + right[r+1] * temp ;
        saved = left[j-r] * temp ;
     }
     N[j] = saved ;
  }

  return;
}




void HB_DersBasisFuns(int p, double start, double incr, double u, int n, double** ders)
{
  int j(p+1), i, r, s1, s2, rk, pk, j1, j2, k;

  double  saved, temp, d;
  std::vector<double>  left(j), right(j);

  vector<vector<double> >  ndu, a;

  //double  ndu[j][j], a[j][j] ;

  ndu.resize(j);
  a.resize(j);
  for(i = 0; i < j; i++)
  {
    ndu[i].resize(j);
    a[i].resize(j);
  }


  ndu[0][0]=1.0 ;
  for(j=1; j<=p; j++)
  {
    left[j] = u-(start-incr*(j-1)) ;
    right[j] = start+incr*j-u ;
    saved = 0.0 ;
    for(r=0; r<j; r++)
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
    s1=0; s2=1; // alternate rows in array a
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




