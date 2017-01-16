#include <iostream>
#include <math.h>
#include "NurbsUtilitiesSURFACE.h"

using namespace std;


/*

  Degree elevate a curve "t" times
  Algorithm A5.9 on pg#206 of the Nurbs Book
  t >= 1

*/

void DegreeElevateSurf1D(NurbsSURFACE* surf1, int t, int dir, NurbsSURFACE* surf2)
{
  if(t<0)
  {
    //Error
    return ;
  }

  if (t == 0)
  {
     surf2->Pw = surf1->Pw;
     surf2->U  = surf1->U;
     surf2->V  = surf1->V;
     surf2->p  = surf1->p;
     surf2->q  = surf1->q;

     return ;
  }

  // create local variables
  CNET Pw, Qw;
  KNOTVECTOR U, V, Uh, Vh;
  DEGREE p, q, ph, qh;

  Pw = surf1->Pw;
  U = surf1->U;
  V = surf1->V;
  p = surf1->p;
  q = surf1->q;

  int mU = U.n-1 ;
  int mV = V.n-1 ;
  int nU = mU-p-1;
  int nV = mV-q-1;

  int i,j,k;

  if (dir == 1)
  {

    Vh = V; // as degree is elevated in U-direction, knot vector(V) and degree(q) in V-direction do not change
    qh = q;


    ph = p+t ; // already defined as type DEGREE in input
    int ph2 = ph/2 ;

    // find distinct knot values present in U
    VectorArray<double> X1;
    findunique(U, X1);

    // find the size of new KNOTVECTOR(Uh) and CPOLYGON(Qw) 
    int mUh1 = mU + t*X1.n;
    int nUh = mUh1-ph-1;

    // resize Qw and Uh to accommodate the increased size
    Uh.setDim(mUh1+1);

    Qw.setDim(nUh+1);
    for(int ii=0;ii<=nUh;ii++)
      Qw[ii].setDim(nV+1);


    //actual computation starts from here
 
    ListArray<VectorArray<double> > bezalfs;
    bezalfs.setDim(ph+1);
    for(i=0;i<(ph+1);i++)
      bezalfs[i].setDim(p+1);

    ListArray<ListArray<CPOINT> > bpts; //bpts(p+1, nV+1)
    bpts.setDim(p+1);
    for(i=0;i<(p+1);i++)
      bpts[i].setDim(nV+1);

    ListArray<ListArray<CPOINT> > ebpts;  //ebpts(ph+1, nV+1)
    ebpts.setDim(ph+1);
    for(i=0;i<(ph+1);i++)
      ebpts[i].setDim(nV+1);

    ListArray<ListArray<CPOINT> > Nextbpts; //Nextbpts(p-1, nV+1)
    Nextbpts.setDim(p-1);
    for(i=0;i<(p-1);i++)
      Nextbpts[i].setDim(nV+1);

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

    int mUh = ph ;
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
  

    for(int col=0;col<=nV;col++)    
      Qw[0][col] = Pw[0][col] ;

    for(i=0; i <= ph ; i++)
      Uh[i] = ua ;


    // Initialize the first Bezier segment
    for(i=0;i<=p ;i++)
      for(int col=0;col<=nV;col++)
        bpts[i][col] = Pw[i][col] ;

    while(b<mU) // Big loop thru knot vector
    {
      i=b ;
      while( b<mU && CompareDoubles(U[b],U[b+1]) )
        b++ ;
      mul = b-i+1 ; 
      mUh += mul+t ;
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
            for(int col=0;col<=nV;col++)
              bpts[k][col] = alphas[k-s]*bpts[k][col] + (1.0-alphas[k-s])*bpts[k-1][col] ;

          for(int col=0;col<=nV;col++)
          Nextbpts[save][col] = bpts[p][col] ;
        }
      } // End of insert knot
    

      for(i=lbz;i<=ph;i++) // Degree elevate Bezier,  only the points lbz,...,ph are used
      {
        for(int col=0;col<=nV;col++)
          ebpts[i][col] = 0.0;
        mpi = min(p,i) ;
        for(j=max(0,i-t); j<=mpi ; j++)
          for(int col=0;col<=nV;col++)
            ebpts[i][col] = ebpts[i][col] + bezalfs[i][j]*bpts[j][col] ;
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
              for(int col=0;col<=nV;col++)
                Qw[i][col] = alf*Qw[i][col] + (1.0-alf)*Qw[i-1][col] ;
	     }
	     if( j>= lbz)
            {
	       if(j-tr <= kind-ph+oldr)
              {
	         gam = (ub-Uh[j-tr])/den ;
                for(int col=0;col<=nV;col++)
                  ebpts[kj][col] = gam*ebpts[kj][col] + (1.0-gam)*ebpts[kj+1][col] ;
	       }
	       else
                for(int col=0;col<=nV;col++)
	         ebpts[kj][col] = bet*ebpts[kj][col] + (1.0-bet)*ebpts[kj+1][col] ;
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
        for(int col=0;col<=nV;col++)
          Qw[cind][col] = ebpts[j][col] ;
        cind++;
      }

      if(b<mU) // Set up for next pass thru loop
      {
        for(j=0;j<r;j++)
          for(int col=0;col<=nV;col++)
            bpts[j][col] = Nextbpts[j][col] ;
        for(j=r;j<=p;j++)
          for(int col=0;col<=nV;col++)
            bpts[j][col] = Pw[b-p+j][col] ;
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
  } // end of if(dir == 1)

  else // dir == 2
  {

    Uh = U; // as degree is elevated in V-direction, knot vector(U) and degree(p) in U-direction do not change
    ph = p;


    qh = q+t ; // already defined as type DEGREE in input
    int qh2 = qh/2 ;

    // find distinct knot values present in V
    VectorArray<double> X1;
    findunique(V, X1);

    // find the size of new KNOTVECTOR(Vh) and CPOLYGON(Qw) 
    int mVh1 = mV + t*X1.n;
    int nVh = mVh1-qh-1;

    // resize Qw and Uh to accommodate the increased size
    Vh.setDim(mVh1+1);

    Qw.setDim(nU+1);
    for(int ii=0;ii<=nU;ii++)
      Qw[ii].setDim(nVh+1);


    //actual computation starts from here
    ListArray<VectorArray<double> > bezalfs; //bezalfs(qh+1,q+1)
    bezalfs.setDim(qh+1);
    for(i=0;i<(qh+1);i++)
      bezalfs[i].setDim(q+1);

    ListArray<ListArray<CPOINT> > bpts; //bpts(nU+1, q+1)
    bpts.setDim(nU+1);
    for(i=0;i<(nU+1);i++)
      bpts[i].setDim(q+1);

    ListArray<ListArray<CPOINT> > ebpts;  //ebpts(nU+1, qh+1)
    ebpts.setDim(nU+1);
    for(i=0;i<(nU+1);i++)
      ebpts[i].setDim(qh+1);

    ListArray<ListArray<CPOINT> > Nextbpts; //Nextbpts(nU+1, q-1)
    Nextbpts.setDim(nU+1);
    for(i=0;i<(nU+1);i++)
      Nextbpts[i].setDim(q-1);

    VectorArray<double> alphas;
    alphas.setDim(q-1);
 

    // Compute Bezier degree elevation coefficients
    double inv, mqi ;
    bezalfs[0][0] = bezalfs[qh][q] = 1.0 ;

    for(i=1;i<=qh2;i++)
    {
      inv= 1.0/Bin(qh,i) ;
      mqi = min(q,i) ;

      for(j=max(0,i-t); j<=mqi; j++)
        bezalfs[i][j] = inv*Bin(q,j)*Bin(t,i-j) ;
    }

    for(i=qh2+1;i<=qh-1 ; i++)
    {
      mqi = min(q,i) ;
      for(j=max(0,i-t); j<=mqi ; j++)
        bezalfs[i][j] = bezalfs[qh-i][q-j] ;
    }

    int mVh = qh ;
    int kind = qh+1 ;
    int r=-1 ; 
    int oldr ;
    int a = q ;
    int b = q+1 ; 
    int cind = 1 ;
    int rbz, lbz ; 
    int mul,save,s;
    double va = V[0] ;
    double vb = 0.0 ;
    double alf;
    int first, last, kj ;
    double den,bet,gam,numer ;
  

    for(int row=0;row<=nU;row++)    
      Qw[row][0] = Pw[row][0] ;

    for(i=0; i <= qh ; i++)
      Vh[i] = va ;


    // Initialize the first Bezier segment
    for(i=0;i<=q ;i++)
      for(int row=0;row<=nU;row++)
        bpts[row][i] = Pw[row][i] ;

    while(b<mV) // Big loop thru knot vector
    {
      i=b ;
      while( b<mV && CompareDoubles(V[b],V[b+1]) )
        b++ ;
      mul = b-i+1 ; 
      mVh += mul+t ;
      vb = V[b] ;
      oldr = r ;
      r = q-mul ;
      // insert V[b] knot "r" times
      if(oldr>0)	lbz = (oldr+2)/2 ;
      else		lbz = 1 ;
      if(r>0)		rbz = qh-(r+1)/2 ;
      else		rbz = qh ;


      if(r>0) // Insert knot to get Bezier segment
      {
        numer = vb-va ;
        for(k=q;k>mul;k--)
          alphas[k-mul-1] = numer/(V[a+k]-va) ;
     
        for(j=1;j<=r;j++)
        {
          save = r-j ; 
          s = mul+j ;
          for(k=q;k>=s;k--)
            for(int row=0;row<=nU;row++)
              bpts[row][k] = alphas[k-s]*bpts[row][k] + (1.0-alphas[k-s])*bpts[row][k-1] ;

          for(int row=0;row<=nU;row++)
          Nextbpts[row][save] = bpts[row][q] ;
        }
      } // End of insert knot
    

      for(i=lbz;i<=qh;i++) // Degree elevate Bezier,  only the points lbz,...,qh are used
      {
        for(int row=0;row<=nU;row++)
          ebpts[row][i] = 0.0;
        mqi = min(q,i) ;
        for(j=max(0,i-t); j<=mqi ; j++)
          for(int row=0;row<=nU;row++)
            ebpts[row][i] = ebpts[row][i] + bezalfs[i][j]*bpts[row][j] ;
      }  
      // End of degree elevating Bezier

      if(oldr>1) // Must remove knot v=V[a] oldr times
      {
        first = kind-2 ; last = kind ;
        den = vb-va ;
        bet = (vb-Vh[kind-1])/den ;
        for(int tr=1; tr<oldr; tr++) // Knot removal loop
        {
	   i = first ; j = last ;
	   kj = j-kind+1 ;
	   while(j-i > tr) // Loop and compute the new control points for one removal step
          {
	     if(i<cind)
            {
	       alf = (vb-Vh[i])/(va-Vh[i]) ;
              for(int row=0;row<=nU;row++)
                Qw[row][i] = alf*Qw[row][i] + (1.0-alf)*Qw[row][i-1] ;
	     }
	     if( j>= lbz)
            {
	       if(j-tr <= kind-qh+oldr)
              {
	         gam = (vb-Vh[j-tr])/den ;
                for(int row=0;row<=nU;row++)
                  ebpts[row][kj] = gam*ebpts[row][kj] + (1.0-gam)*ebpts[row][kj+1] ;
	       }
	       else
                for(int row=0;row<=nU;row++)
	         ebpts[row][kj] = bet*ebpts[row][kj] + (1.0-bet)*ebpts[row][kj+1] ;
	     }
	     i++ ; j--; kj-- ;
	   }
	   first-- ; last++ ;
         }
       }
      // end of removing knot v=V[a]
      if(a!=q) // load the knot va
      {
        for(i=0;i<qh-oldr; i++)
        {
          Vh[kind] = va ;
          kind++;
        }
      }
      for(j=lbz; j<=rbz ; j++) // load control points into Qw
      {
        for(int row=0;row<=nU;row++)
          Qw[row][cind] = ebpts[row][j] ;
        cind++;
      }

      if(b<mV) // Set up for next pass thru loop
      {
        for(j=0;j<r;j++)
          for(int row=0;row<=nU;row++)
            bpts[row][j] = Nextbpts[row][j] ;
        for(j=r;j<=q;j++)
          for(int row=0;row<=nU;row++)
            bpts[row][j] = Pw[row][b-q+j] ;
        a = b ; 
        b++ ;
        va = vb ;
      }
      else
      {
        for(i=0;i<=qh;i++)
         Vh[kind+i] = vb ;
      }
    }
  } // end of else //if(dir == 1) loop

    // assign values to member variables of surf2
    surf2->Pw = Qw;
    surf2->U = Uh;
    surf2->V = Vh;
    surf2->p = ph;
    surf2->q = qh;

}
