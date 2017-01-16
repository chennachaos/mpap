#include <iostream>
#include <math.h>
#include "NurbsUtilitiesSURFACE.h"

using namespace std;


/*
  Refine the Surface knot vector in one direction
  Algorithm A5.5 on pg#167 of the Nurbs Book
  
  X the new knots to insert in the knot vector
*/

void RefineSurfKnotVector1D(NurbsSURFACE* surf1, KNOTVECTOR& X, int dir, NurbsSURFACE* surf2)
{
  // create local variables
  CNET Pw, Qw;
  KNOTVECTOR U, V, Ubar, Vbar;
  DEGREE p, q;

  Pw = surf1->Pw;
  U = surf1->U;
  V = surf1->V;
  p = surf1->p;
  q = surf1->q;

  int mU = U.n-1 ;
  int mV = V.n-1 ;
  int nU = mU-p-1;
  int nV = mV-q-1;

  int r = X.n-1 ;
  int i, j, k, l, ii, jj, col, row, ind ;
  
  double alpha;

 
  if( dir == 1)
  {
    // resize Qw and Ubar to accommodate the increased size
    int nUbar = nU + r + 1;
    int mUbar = mU + r + 1;
    
    Qw.setDim(nUbar+1);
    for(ii=0;ii<=nUbar;ii++)
      Qw[ii].setDim(nV+1);

    Ubar.setDim(mUbar+1);
   
    int a = FindSpan(&(U[0]), U.n, p, X[0]) ;
    int b = FindSpan(&(U[0]), U.n, p, X[r]) ;
    b=b+1;
    
    Vbar = V;  
    
    //Save unaltered control points
    for(col=0;col<=nV;col++)
    {
      for(k=0; k<=a-p ; k++)
        Qw[k][col] = Pw[k][col] ;
      for(k=b-1; k<=nU; k++)
        Qw[k+r+1][col] = Pw[k][col] ;
    }
  
    for(j=0; j<=a ; j++)
      Ubar[j] = U[j] ;
    for(j=b+p; j<=mU; j++)
      Ubar[j+r+1] = U[j] ;

    i = b+p-1 ; 
    k = b+p+r ;

    for(j=r; j>=0 ; j--)
    {
      while(X[j] <= U[i] && i>a)
      {
        for(col=0;col<=nV;col++)
          Qw[k-p-1][col] = Pw[i-p-1][col] ;

        Ubar[k] = U[i] ;
        k = k-1;
        i = i-1;
      }

      for(col=0;col<=nV;col++)    
        Qw[k-p-1][col] = Qw[k-p][col] ;
    
      for(jj=1; jj<=p ; jj++)
      {
        ind = k-p+jj ;
        alpha = Ubar[k+jj] - X[j] ;
        if(CompareDoubles(fabs(alpha), 0.0))
        {
          for(col=0;col<=nV;col++)
            Qw[ind-1][col] = Qw[ind][col] ;
        }
        else
        {
	   alpha = alpha/(Ubar[k+jj]-U[i-p+jj]) ;
          for(col=0;col<=nV;col++)
            Qw[ind-1][col] = alpha*Qw[ind-1][col] + (1.0-alpha)*Qw[ind][col] ;
        }
      }
      Ubar[k] = X[j] ;
      k = k-1;
    }
    // assign values to member variables of surf2
    surf2->Pw = Qw;
    surf2->U = Ubar;
    surf2->V = V;
    surf2->p = p;
    surf2->q = q;

  } // if(dir == 1) ends here

  
  else // if dir != 1 then dir == 2
  {

    // resize Qw and Vbar to accommodate the increased size
    int nVbar = nV + r + 1;
    int mVbar = mV + r + 1;

    Vbar.setDim(mVbar+1); Vbar.zero();
    Qw.setDim(nU+1);
    for(ii=0;ii<=nU;ii++)
      Qw[ii].setDim(nVbar+1);

    int a = FindSpan(&(V[0]), V.n, q, X[0]) ;
    int b = FindSpan(&(V[0]), V.n, q, X[r]) ;
    b=b+1;
    
    Ubar = U;  
    
    //Save unaltered control points
    for(row=0;row<=nU;row++)
    {
      for(k=0; k<=a-q ; k++)
        Qw[row][k] = Pw[row][k] ;
      for(k=b-1; k<=nV; k++)
        Qw[row][k+r+1] = Pw[row][k] ;
    }
  
    for(j=0; j<=a ; j++)
      Vbar[j] = V[j] ;
    for(j=b+q; j<=mV; j++)
      Vbar[j+r+1] = V[j] ;

    i = b+q-1 ;
    k = b+q+r ;

    for(j=r; j>=0 ; j--)
    {
      while(X[j] <= V[i] && i>a)
      {
        for(row=0;row<=nU;row++)
          Qw[row][k-q-1] = Pw[row][i-q-1] ;

        Vbar[k] = V[i] ;
        k = k-1;
        i = i-1;
      }

      for(row=0;row<=nU;row++)    
        Qw[row][k-q-1] = Qw[row][k-q] ;
    
      for(l=1; l<=q ; l++)
      {
        ind = k-q+l ;
        alpha = Vbar[k+l] - X[j] ;
        if(CompareDoubles(fabs(alpha), 0.0))
        {
          for(row=0;row<=nU;row++)
            Qw[row][ind-1] = Qw[row][ind] ;
        }
        else
        {
	   alpha = alpha/(Vbar[k+l]-V[i-q+l]) ;
          for(row=0;row<=nU;row++)
            Qw[row][ind-1] = alpha*Qw[row][ind-1] + (1.0-alpha)*Qw[row][ind] ;
        }
      }
      Vbar[k] = X[j] ;
      k = k-1;
    }

    // assign values to member variables of surf2
    surf2->Pw = Qw;
    surf2->U = U;
    surf2->V = Vbar;
    surf2->p = p;
    surf2->q = q;
  } // else ends here
}
