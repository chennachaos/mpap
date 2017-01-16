#include <iostream>
#include <math.h>
#include "NurbsUtilitiesCURVE.h"

using namespace std;



/*

  Degree elevate a curve "t" times
  Algorithm A5.9 on pg#206 of the Nurbs Book
  t >= 1

*/

void DegreeElevateCurve(NurbsCURVE* curv1, int t, NurbsCURVE* curv2)
{
  if(t<0)
  {
    //Error
    return ;
  }

  if (t == 0)
  {
    curv1 = curv2;
    return ;
  }

  CPOLYGON Pw, Qw;
  KNOTVECTOR U, Uh;
  DEGREE p, ph;
  Pw = curv1->Pw;
  U = curv1->U;
  p = curv1->p;

  int i,j,k ;
  int n = Pw.n-1 ;
  int m = U.n-1 ;
  
  ph = p+t ; // already defined as DEGREE in input
  int ph2 = ph/2 ;

  // find distinct knot values present in U
  VectorArray<double> X1;
  findunique(U, X1);

  // find the size of new KNOTVECTOR(Uh) and CPOLYGON(Qw) 
  int mh1 = m + t*X1.n;
  int nh = mh1-ph-1;

  // resize Qw and Uh to accommodate the increased size
  Qw.setDim(nh+1);
  Uh.setDim(mh1+1);


  ListArray<VectorArray<double> > bezalfs;
  bezalfs.setDim(ph+1);
  for(i=0;i<(ph+1);i++)
    bezalfs[i].setDim(p+1);

  ListArray<CPOINT> bpts;
    bpts.setDim(p+1);
  ListArray<CPOINT> ebpts;
    ebpts.setDim(ph+1);
  ListArray<CPOINT> Nextbpts;
    Nextbpts.setDim(p-1);
  VectorArray<double> alphas;
    alphas.setDim(p-1);


  // Compute Bezier degree elevation coefficients
  double inv, mpi ;
  bezalfs[0][0] = bezalfs[ph][p] = 1.0 ;

  for(i=1;i<=ph2;i++)
  {
    inv= 1.0/Bin(ph,i) ;
    mpi = min(p,i) ;

    for(j=max(0,i-t); j<=mpi; j++)
      bezalfs[i][j] = inv*Bin(p,j)*Bin(t,i-j) ;
  }

  for(i=ph2+1;i<=ph-1 ; i++)
  {
    mpi = min(p,i) ;
    for(j=max(0,i-t); j<=mpi ; j++)
      bezalfs[i][j] = bezalfs[ph-i][p-j] ;
  }

  int mh = ph ;
  int kind = ph+1 ;
  int r=-1 ; 
  int oldr ;
  int a = p ;
  int b = p+1 ; 
  int cind = 1 ;
  int rbz, lbz ; 
  int mul,save,s;
  double ua = U[0] ;
  double ub = 0.0 ;
  double alf;
  int first, last, kj ;
  double den,bet,gam,numer ;
  

  Qw[0] = Pw[0] ;
  for(i=0; i <= ph ; i++)
    Uh[i] = ua ;

  // Initialize the first Bezier segment
  for(i=0;i<=p ;i++) 
    bpts[i] = Pw[i] ;

  while(b<m) // Big loop thru knot vector
  {
    i=b ;
    while( b<m && CompareDoubles(U[b],U[b+1]) )
      b++ ;
    mul = b-i+1 ; 
    mh += mul+t ;
    ub = U[b] ;
    oldr = r ;
    r = p-mul ;
    // insert U[b] knot "r" times
    if(oldr>0)	lbz = (oldr+2)/2 ;
    else		lbz = 1 ;
    if(r>0)		rbz = ph-(r+1)/2 ;
    else		rbz = ph ;

    if(r>0) // Insert knot to get Bezier segment
    {
      numer = ub-ua ;
      for(k=p;k>mul;k--)
        alphas[k-mul-1] = numer/(U[a+k]-ua) ;
     
      for(j=1;j<=r;j++)
      {
        save = r-j ; 
        s = mul+j ;
	 for(k=p;k>=s;k--)
          bpts[k] = alphas[k-s]*bpts[k] + (1.0-alphas[k-s])*bpts[k-1] ;
        Nextbpts[save] = bpts[p] ;
      }
    } // End of insert knot
    

    for(i=lbz;i<=ph;i++) // Degree elevate Bezier,  only the points lbz,...,ph are used
    {
      ebpts[i] = 0.0;
      mpi = min(p,i) ;
      for(j=max(0,i-t); j<=mpi ; j++)
        ebpts[i] = ebpts[i] + bezalfs[i][j]*bpts[j] ;
    }  
    // End of degree elevating Bezier

    if(oldr>1) // Must remove knot u=U[a] oldr times
    {
      first = kind-2 ; last = kind ;
      den = ub-ua ;
      bet = (ub-Uh[kind-1])/den ;
      for(int tr=1; tr<oldr; tr++) // Knot removal loop
      {
	 i = first ; j = last ;
	 kj = j-kind+1 ;
	 while(j-i > tr) // Loop and compute the new control points for one removal step
        {
	   if(i<cind)
          {
	     alf = (ub-Uh[i])/(ua-Uh[i]) ;
	     Qw[i] = alf*Qw[i] + (1.0-alf)*Qw[i-1] ;
	   }
	   if( j>= lbz)
          {
	     if(j-tr <= kind-ph+oldr)
            {
	       gam = (ub-Uh[j-tr])/den ;
	       ebpts[kj] = gam*ebpts[kj] + (1.0-gam)*ebpts[kj+1] ;
	     }
	     else
	       ebpts[kj] = bet*ebpts[kj] + (1.0-bet)*ebpts[kj+1] ;
	   }
	   i++ ; j--; kj-- ;
	 }
	 first-- ; last++ ;
       }
     }
    // end of removing knot u=U[a]
    if(a!=p) // load the knot ua
    {
      for(i=0;i<ph-oldr; i++)
      {
        Uh[kind] = ua ;
        kind++;
      }
    }
    for(j=lbz; j<=rbz ; j++) // load control points into Qw
    {
      Qw[cind] = ebpts[j] ;
      cind++;
    }

    if(b<m) // Set up for next pass thru loop
    {
      for(j=0;j<r;j++)
	 bpts[j] = Nextbpts[j] ;
      for(j=r;j<=p;j++)
	 bpts[j] = Pw[b-p+j] ;
      a = b ; 
      b++ ;
      ua = ub ;
    }
    else
    {
      for(i=0;i<=ph;i++)
	 Uh[kind+i] = ub ;
    }
  }

  curv2->Pw = Qw;
  curv2->U = Uh;
  curv2->p = ph;

    return ;
}
