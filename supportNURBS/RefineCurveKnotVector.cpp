#include <iostream>
#include <math.h>
#include "NurbsUtilitiesCURVE.h"

using namespace std;


/*
  Refine the curve knot vector
  Algorithm A5.4 on pg#165 of the Nurbs Book
  
  X the new knots to insert in the knot vector
*/

void RefineCurveKnotVector(NurbsCURVE* curv1, KNOTVECTOR& X, NurbsCURVE* curv2)
{
  CPOLYGON Pw, Qw;
  KNOTVECTOR U, Ubar;
  DEGREE p;
  Pw = curv1->Pw;
  U = curv1->U;
  p = curv1->p;

  double alpha;

  int mU = U.n-1 ;
  int nU = mU-p-1;

  int  r = X.n-1 ;
  int  a, b, j, l, ind ;

  int nUbar = nU + r+1;
  int mUbar = mU + r+1;
  // resize Qw and Ubar to accommodate the increased size
  Qw.setDim(nUbar+1);
  Ubar.setDim(mUbar+1);


  a = FindSpan(&(U[0]), U.n, p, X[0]) ;
  b = FindSpan(&(U[0]), U.n, p, X[r]) ;
  b=b+1;

  for(j=0; j<=a-p ; j++)
    Qw[j] = Pw[j] ;
  for(j=b-1; j<=nU; j++)
    Qw[j+r+1] = Pw[j] ;
  for(j=0; j<=a ; j++)
    Ubar[j] = U[j] ;
  for(j=b+p; j<=mU; j++)
    Ubar[j+r+1] = U[j] ;

  int i = b+p-1 ; 
  int k = b+p+r ;

  for(j=r; j>=0 ; j--)
  {
    while(X[j] <= U[i] && i>a)
    {
      Qw[k-p-1] = Pw[i-p-1] ;
      Ubar[k] = U[i] ;
      k = k-1;
      i = i-1;
    }
    Qw[k-p-1] = Qw[k-p] ;
    
    for(l=1; l<=p ; l++)
    {
      ind = k-p+l ;
      alpha = Ubar[k+l] - X[j] ;

      if(CompareDoubles(fabs(alpha), 0.0))
	 Qw[ind-1] = Qw[ind] ;
      else
      {
	 alpha = alpha/(Ubar[k+l]-U[i-p+l]) ;
        Qw[ind-1] = alpha*Qw[ind-1] + (1.0-alpha)*Qw[ind] ;
      }
    }
    Ubar[k] = X[j] ;
    k = k-1;
  }

  curv2->Pw = Qw;
  curv2->U = Ubar;
  curv2->p = p;

    return ;
}
