#include <iostream>
#include <math.h>
#include "NurbsUtilitiesCURVE.h"

using namespace std;



int RepamCurveNelem(NurbsCURVE* curv1, int p, int Nel1, NurbsCURVE* curv2)
{
    int mp = 2*p + Nel1;
    int np = mp-p-1, ii, jj, ll;

    double invNel1, invnp, temp, temp1, xlen;
    invNel1 = 1.0/Nel1;
    invnp   = 1.0/np;


    KNOTVECTOR U, UU;
    U.setDim(mp+1);
    UU.setDim(np+1);
    U.zero();    UU.zero();

    EPOINT EP1, EP2;
    EP1 = curv1->Pw[0].CalcEuclid();
    EP2 = curv1->Pw[curv1->Pw.n-1].CalcEuclid();

    xlen = CalcDist(EP1, EP2);

   //  cout << " xlen & ylen " << xlen << '\t' << ylen << endl;

/*
    // find knot vector in 1st direction
    for(int ii=0;ii<(Nel1-1);ii++)
       U[p+1+ii] = invNel1 * (ii+1);
    for(int ii=0;ii<=p;ii++)
       U[mp-ii] = 1.0;
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


/*
    cout << '\t' << " np " << np << endl;
    cout << endl;
    cout << '\t' << U << endl;
    cout << endl;
*/
    ListArray<CPOINT> Qw;
    Qw.setDim(np+1);

    temp1 = invnp*xlen;

    for(ii=0;ii<=np;ii++)
    {
       Qw[ii]   = 0.0;
       Qw[ii].x = temp1 * ii;
       Qw[ii].w = 1.0;
    }

      curv2->Pw.setDim(Qw.n);
      curv2->U.setDim(U.n);

      curv2->Pw = Qw;
      curv2->U  = U;
      curv2->p  = p;

      curv2->initializeBCdata();


      double u0=0.3, val11;

      VectorArray<double> ubar1;

      int nn=np-p;
      ubar1.setDim(nn);
      EP1 = 0.0;
      val11 = xlen/(nn+1);

      // find the parameter value
      for(ii=0;ii<nn;ii++)
      {
         EP1.x = (ii+1) * val11;
         ubar1[ii] = curv1->CurvePointInverse(EP1, u0);
      }

      for(ii=0;ii<nn;ii++)
        curv2->U[ii+p+1] = ubar1[ii];

/*
    cout << endl;
    cout << '\t' << curv2->U << endl;
    cout << endl;
*/

 return 0;

}
