#include <iostream>
#include <math.h>
#include "NURBS_1D.h"

using namespace std;


double OneBasisAlg1(double* U, int nU, DEGREE p, int i, double u)
{
	//cout << " i and deg " << '\t' << i << '\t' << p << endl;

    if( (CompareDoubles(u, U[0]) && (i == 0) ) || ( CompareDoubles(u, U[nU-1])  && (i == nU-p-2)) )
      return 1.0;

    if( u<U[i] || u>=U[i+p+1] )
      return 0.0;


    VectorArray<double> N; N.setDim(p+1);
    double saved, Uleft, Uright, temp;
    int j, k;

    for(j=0;j<=p;j++)
    {
        if(u>=U[i+j] && u<U[i+j+1])
            N[j] = 1.0;
        else
            N[j] = 0.0;
    }

    for(k=1;k<=p;k++)
    {
        if( CompareDoubles(N[0], 0.0) )
           saved = 0.0;
        else
           saved = ((u-U[i])*N[0])/(U[i+k]-U[i]);

        for(j=0;j<p-k+1;j++)
        {
            Uleft  = U[i+j+1];
            Uright = U[i+j+k+1];
            if( CompareDoubles(N[j+1], 0.0) )
            {
                N[j] = saved;
                saved = 0.0;
            }
            else
            {
                temp = N[j+1]/(Uright-Uleft);
                N[j] = saved + (Uright-u)*temp;
                saved = (u-Uleft)*temp;
            }
        }
    }

    return N[0];
}

