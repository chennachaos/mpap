#include <iostream>
#include <math.h>
#include "NurbsUtilitiesSURFACE.h"

using namespace std;


// find knot vector using reparameterization

void FindKnotVector(NurbsSURFACE* surf1, int dir, int Nsub, KNOTVECTOR& XX)
{
  VectorArray<double> X1;

  NurbsCURVE curv1, curv2;

  EPOINT EP1, EP2, EP3, EP4, EP5, EP6, EP7, EP8;

  CPOINT CP1, CP2, CP3, CP4, CP5, CP6, CP7, CP8;

  double ubar1=0.0, ubar2=0.0, incr=0.0, xx=0.0, yy=0.0, dummy=0.0, val;
  double dist1=0.0, dist2=0.0, u0=0.0, temp=0.0, norm1, norm2, uw;

  ListArray<EPOINT> derEP1, C1, C2;

  int N = 25, count = 0, span, size;

   ListArray<VectorArray<double> > ders1;

   VectorArray<double> dN_du;


  if(dir == 1) 
  {
     findunique(surf1->U, X1);

     // get the v=0.0 isocurve
     surf1->CurveOnSurf(2, 0.0, &curv1);

     // get the v=1.0 isocurve
     surf1->CurveOnSurf(2, 1.0, &curv2);


     dN_du.setDim(surf1->p+1);	dN_du.zero();
  }
  else 
  {
     findunique(surf1->V, X1);

     // get the u=0.0 isocurve
     surf1->CurveOnSurf(1, 0.0, &curv1);

     // get the u=1.0 isocurve
     surf1->CurveOnSurf(1, 1.0, &curv2);

     dN_du.setDim(surf1->q+1);	dN_du.zero();
  }


     val =  pow(2.0,Nsub);
     size = (int) val - 1;

     XX.setDim((X1.n-1) * size);

     C1.setDim(N-2); C2.setDim(N-2);

     for(int jj=0;jj<X1.n-1;jj++)
     {
           // calculate the end points at the knot values
           EP1 = curv1.CurvePoint(X1[jj]).CalcEuclid();
           EP2 = curv1.CurvePoint(X1[jj+1]).CalcEuclid();

           // compute the paramter value for the second curve
           EP3 = curv2.CurvePoint(X1[jj]).CalcEuclid();
           EP4 = curv2.CurvePoint(X1[jj+1]).CalcEuclid();

    /*
           EP1.print2screen();
           cout << endl;
           EP2.print2screen();
           cout << endl;
           EP3.print2screen();
           cout << endl;
           EP4.print2screen();
           cout << endl;
    */
           incr = (X1[jj+1] - X1[jj])/N;

           for(int ll=1;ll<(N-1);ll++)
           {
              temp = X1[jj]+ll*incr;
              C1[ll-1] = curv1.CurvePoint(temp).CalcEuclid();
              C2[ll-1] = curv2.CurvePoint(temp).CalcEuclid();
           }

        for(int ii=0;ii<size;ii++)
        {
           // compute the paramter value for the first curve
           dummy = ((ii+1)/val);

           EP5 = EP1 + dummy * (EP2 - EP1);

           // EP5.print2screen();
           // cout << endl;

           // find the initial guess
           dist1 = CalcDist(C1[0], EP5);
           u0 = X1[jj]+incr;

           dist2 = dist1;
           for(int ll=1;ll<N-1;ll++)
           {
              //C1[ll-1].print2screen();
              temp = X1[jj]+(ll+1)*incr;
              dist2 = CalcDist(C1[ll-1], EP5);

              if(doubleLess(dist2, dist1) )
              {
                 dist1 = dist2;
                 u0 = temp;
              }
           }

           //  cout << "         u0 : " << u0 << endl;

           // find the parameter value
           ubar1 = curv1.CurvePointInverse(EP5, u0);

       /*
           cout << endl;
           cout << "         ubar1 : " << ubar1 << endl;
           cout << endl;
       */

           curv1.CurveDerPointRat(ubar1, 1, derEP1);
           norm1 = derEP1[1].Norm();


           EP7 = EP3 + dummy * (EP4 - EP3);

           // EP7.print2screen();

           dist1 = CalcDist(C2[0], EP7);
           u0 = X1[jj]+incr;

           dist2 = dist1;
           for(int ll=1;ll<N-1;ll++)
           {
              //C2[ll-1].print2screen();
              temp = X1[jj]+(ll+1)*incr;
              dist2 = CalcDist(C2[ll-1], EP7);

              if(doubleLess(dist2, dist1) )
              {
                dist1 = dist2;
                u0 = temp;
              }
           }

           ubar2 = curv2.CurvePointInverse(EP7, u0);

       /*
           cout << endl;
           cout << "         ubar2 : " << ubar2 << endl;
           cout << endl;
       */
           curv2.CurveDerPointRat(ubar2, 1, derEP1);
           norm2 = derEP1[1].Norm();

           // compute and store the weighted paramter value

           uw = (norm1 * ubar1 + norm2 * ubar2)/(norm1+norm2);
           XX[count++] = uw;
          // XX[count++] = (ubar1 + ubar2)/2;

       /*
           cout << endl;
           cout << "         uw : " << uw << endl;
           cout << endl;
       */

        }
     }
       /*
          for(int ii=0;ii<XX.n;ii++)
            cout << '\t' << XX[ii] << endl;
          cout << endl;
       */

    SortArray(XX);

  return;

}
