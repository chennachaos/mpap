#include <iostream>
#include <math.h>
#include "NurbsUtilitiesSURFACE.h"

using namespace std;



int RepamSurf2DNelem(NurbsSURFACE* surf0, int p, int q, int Nel1, int Nel2, NurbsSURFACE* surf1)
{
    int ii, jj, nn, ll;

    int mp = 2*p + Nel1;
    int mq = 2*q + Nel2;
    int np = mp-p-1;
    int nq = mq-q-1;

    double invNel1, invNel2, invnp, invnq, temp, temp1, temp2, xlen, ylen;
    invNel1 = 1.0/Nel1;
    invNel2 = 1.0/Nel2;
    invnp   = 1.0/np;
    invnq   = 1.0/nq;


    KNOTVECTOR U, V, UU, VV;
    U.setDim(mp+1);
    V.setDim(mq+1);
    UU.setDim(np+1);
    VV.setDim(nq+1);
    U.zero();     V.zero();
    UU.zero();    VV.zero();

    EPOINT EP1, EP2;
    EP1 = surf0->Pw[0][0].CalcEuclid();
    EP2 = surf0->Pw[surf0->Pw.n-1][0].CalcEuclid();

    xlen = CalcDist(EP1, EP2);

    EP2 = surf0->Pw[0][surf0->Pw[0].n-1].CalcEuclid();

    ylen = CalcDist(EP1, EP2);

   //  cout << " xlen & ylen " << xlen << '\t' << ylen << endl;

/*
    // find knot vector in 1st direction
    for(int ii=0;ii<(Nel1-1);ii++)
       U[p+1+ii] = invNel1 * (ii+1);
    for(int ii=0;ii<=p;ii++)
       U[mp-ii] = 1.0;

    // find knot vector in 2nd direction
    for(int ii=0;ii<(Nel2-1);ii++)
       V[q+1+ii] = invNel2 * (ii+1);
    for(int ii=0;ii<=q;ii++)
       V[mq-ii] = 1.0;
*/


    // find knot vector in 1st direction
    for(ii=0;ii<=np;ii++)
       UU[ii] = invnp*ii;
    for(ii=0;ii<=p;ii++)
       U[mp-ii] = 1.0;
    for(jj=1;jj<=(np-p);jj++)
    {
       for(ii=jj;ii<=(jj+p-1);ii++)
         U[jj+p] += UU[ii];
       U[jj+p] = U[jj+p]/p;
    }

    // find knot vector in 2nd direction
    for(ii=0;ii<=nq;ii++)
       VV[ii] = invnq*ii;
    for(ii=0;ii<=q;ii++)
       V[mq-ii] = 1.0;
    for(jj=1;jj<=(nq-q);jj++)
    {
       for(ii=jj;ii<=(jj+q-1);ii++)
         V[jj+q] += VV[ii];
       V[jj+q] = V[jj+q]/q;
    }



//
    cout << '\t' << " np & nq " << np << '\t' << nq << endl;
    cout << endl;
    cout << '\t' << U << endl;
    cout << endl;
    cout << '\t' << V << endl;
    cout << endl;
//
    ListArray<ListArray<CPOINT> > Qw;
    Qw.setDim(np+1);
    for(ii=0;ii<=np;ii++)
      Qw[ii].setDim(nq+1);

    temp1 = invnp*xlen;
    temp2 = invnq*ylen;

    for(ii=0;ii<=np;ii++)
    {
       temp = temp1 * ii;
       for(jj=0;jj<=nq;jj++)
       {
         Qw[ii][jj]   = 0.0;
         Qw[ii][jj].x = temp;
         Qw[ii][jj].y = temp2 * jj;
         Qw[ii][jj].w = 1.0;
       }
    }

      surf1->Pw.setDim(Qw.n);
      for(ll=0;ll<Qw[0].n;ll++)
        surf1->Pw[ll].setDim(Qw[0].n);

      surf1->U.setDim(U.n);
      surf1->V.setDim(V.n);

      surf1->Pw = Qw;
      surf1->U  = U;
      surf1->V  = V;
      surf1->p  = p;
      surf1->q  = q;

      surf1->initializeBCdata();

   if(p > 1)
   {
      // DIRECTION #1
      // ----------------
      double u0=0.3, val11;
      NurbsCURVE curv1;
      VectorArray<double> ubar1;

      // get the v=0.0 isocurve
      surf1->CurveOnSurf(2, 0.0, &curv1);

      nn=np-p;
      ubar1.setDim(nn);
      EP1 = 0.0;
      val11 = xlen/(nn+1);

      // find the parameter value
      for(ii=0;ii<nn;ii++)
      {
         EP1.x = (ii+1) * val11;
         ubar1[ii] = curv1.CurvePointInverse(EP1, u0);
      }

      for(ii=0;ii<nn;ii++)
        surf1->U[ii+p+1] = ubar1[ii];

      cout << endl;
      cout << '\t' << ubar1 << endl;
   }

   if(q > 1)
   {
      // DIRECTION #2
      // ----------------
      double u0=0.3, val11;
      NurbsCURVE curv1;
      VectorArray<double> ubar1;

      // get the u=0.0 isocurve
      surf1->CurveOnSurf(1, 0.0, &curv1);

      nn=nq-q;
      ubar1.setDim(nn);
      val11 = ylen/(nn+1);

      // find the parameter value
      for(ii=0;ii<nn;ii++)
      {
         EP1.y = (ii+1) * val11;
         ubar1[ii] = curv1.CurvePointInverse(EP1, u0);
      }

      for(ii=0;ii<nn;ii++)
        surf1->V[ii+q+1] = ubar1[ii];

      cout << endl;
      cout << '\t' << ubar1 << endl;
   }

    cout << endl;
    cout << '\t' << surf1->U << endl;
    cout << endl;
    cout << '\t' << surf1->V << endl;
    cout << endl;



 return 0;

}
