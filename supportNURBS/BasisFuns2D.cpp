#include <iostream>
#include <math.h>
#include "NURBS_1D.h"

using namespace std;



void BasisFuns2D(double* U, int Un, int p, double* V, int Vn, int q, double u, double v, double* NN)
{
  vector<double> N(p+1), M(q+1);

  BasisFuns(U, Un, p, u, &N[0]);

  BasisFuns(V, Vn, q, v, &M[0]);

  int loc_num = 0, ii, jj;

  for(jj=0;jj<=q;jj++)
  {
    for(ii=0;ii<=p;ii++)
    {
      NN[loc_num] = N[ii] * M[jj];
      loc_num++;
    }
  }

  return;  
}
