#include <iostream>
#include <math.h>
#include "NurbsUtilitiesCURVE.h"

using namespace std;


/*
  Insert a knot ubar, r times
  Algorithm A5.1 on pg#151 of the Nurbs Book
*/

void CurveKnotIns(NurbsCURVE* curv1, double ubar, int r, NurbsCURVE* curv2)
{
  CPOLYGON Pw, Qw;
  KNOTVECTOR UP, UQ;
  DEGREE p;
  Pw = curv1->Pw;
  UP = curv1->U;
  p = curv1->p;

  int np = Pw.n-1 ;
  int mp = UP.n-1 ;
  int nq = np+r;
  int mq = mp+r;
  double alpha;
  int L, ii, jj, j;
  ListArray<CPOINT> Rw;
  Rw.setDim(p+1);

  // resize Qw and UQ to accommodate the increased size
  Qw.setDim(nq+1);
  UQ.setDim(mq+1);


  int s = FindMult(&(UP[0]), UP.n, p, ubar); // multiplicity of the knot "ubar" in UP

  int k = FindSpan(&(UP[0]), UP.n, p, ubar);

  //Load New Vector
  for(ii=0;ii<=k;ii++)
    UQ[ii] = UP[ii];
  
  for(ii=1;ii<=r;ii++)
    UQ[k+ii] = ubar;
  
  for(ii=k+1;ii<=mp;ii++)
    UQ[ii+r] = UP[ii];

  //Save unaltered control points
  for(ii=0;ii<=k-p;ii++)
    Qw[ii]=Pw[ii];

  for(ii=k-s;ii<=np;ii++)
    Qw[ii+r] = Pw[ii];

  for(ii=0;ii<=p-s;ii++)
    Rw[ii] = Pw[k-p+ii];

  // Insert knot 'r' times
  for(j=1;j<=r;j++)
  {
    L = k-p+j;
    for(ii=0;ii<=p-j-s;ii++)
    {
      alpha = (ubar-UP[L+ii])/(UP[ii+k+1]-UP[L+ii]);
      Rw[ii] = alpha*Rw[ii+1] + (1.0-alpha)*Rw[ii];
    }
    Qw[L] = Rw[0];
    Qw[k+r-j-s] = Rw[p-j-s];
  }
  // Load remaining control points
  for(ii=L+1;ii<k-s;ii++)
    Qw[ii]=Rw[ii-L];

  curv2->Pw = Qw;
  curv2->U = UQ;
  curv2->p = p;

}
