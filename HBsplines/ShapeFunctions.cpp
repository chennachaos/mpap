
#include "Debug.h"
#include "ShapeFunctions.h"



void ShapeFunctions::initialise(int p)
{
    int ind = p+1;

    N.resize(ind);          N.setZero();
    dN_dx.resize(ind);      dN_dx.setZero();
    d2N_dx2.resize(ind);    d2N_dx2.setZero();

    return;
}


void ShapeFunctions::initialise(int p, int q)
{
    int ind = (p+1)*(q+1); 

    N.resize(ind);          N.setZero();
    dN_dx.resize(ind);      dN_dx.setZero();
    d2N_dx2.resize(ind);    d2N_dx2.setZero();
    dN_dy.resize(ind);      dN_dy.setZero();
    d2N_dy2.resize(ind);    d2N_dy2.setZero();
    
    return;
}


void ShapeFunctions::initialise(int p, int q, int r)
{
    int ind = (p+1)*(q+1)*(r+1);

    N.resize(ind);          N.setZero();

    dN_dx = N;
    dN_dy = N;
    dN_dz = N;

    d2N_dx2 = N;
    d2N_dy2 = N;
    d2N_dz2 = N;

    //dN_dx.resize(ind);      dN_dx.setZero();
    //d2N_dx2.resize(ind);    d2N_dx2.setZero();
    //dN_dy.resize(ind);      dN_dy.setZero();
    //d2N_dy2.resize(ind);    d2N_dy2.setZero();
    //dN_dz.resize(ind);      dN_dz.setZero();
    //d2N_dz2.resize(ind);    d2N_dz2.setZero();

    return;
}

void ShapeFunctions::zero()
{
  N.setZero();
  dN_dx.setZero();   d2N_dx2.setZero();
  dN_dy.setZero();   d2N_dy2.setZero();
  dN_dz.setZero();   d2N_dz2.setZero();

  return;
}



void ShapeFunctions::Compute()
{
/*
  int ROWS = 3, ii, jj, count;

    double** ders1 = new double*[ROWS];
    double** ders2 = new double*[ROWS];

    for(ii=0;ii<ROWS;ii++)
    {
       ders1[ii] = new double[degree[0]+1];
       ders2[ii] = new double[degree[1]+1];
    }

    double  fact1, fact2, fact3, fact4, fact5, fact6, uu, vv;
    
    fact1 = Jacobian[0]*Jacobian[0];
    fact2 = Jacobian[1]*Jacobian[1];

    for(gp2=0;gp2<HbsData->getNGP(1);gp2++)
    {
       vv = val3 * gausspoints2[gp2] + val4;
       
       for(gp1=0;gp1<HbsData->getNGP(0);gp1++)
       {
          uu = val1 * gausspoints1[gp1] + val2;

          HB_DersBasisFuns(degree[0], knots[0][0], knots[0][2], uu, 2, ders1);
          HB_DersBasisFuns(degree[1], knots[1][0], knots[1][2], vv, 2, ders2);

          count = 0;
          for(jj=0; jj<=degree[1]; jj++)
          {
             fact3 = ders2[0][jj]/Jacobian[0];
             fact4 = ders2[1][jj]/Jacobian[1];

             fact5 = ders2[0][jj]/fact1;
             fact6 = ders2[2][jj]/fact2;
       
             for(ii=0; ii<=degree[0]; ii++)
             {
                N[count]       =  ders1[0][ii] * ders2[0][jj];
       
                dN_dx[count]   =  ders1[1][ii] * fact3;
                dN_dy[count]   =  ders1[0][ii] * fact4;

                d2N_dx2[count] =  ders1[2][ii] * fact5;
                d2N_dy2[count] =  ders1[0][ii] * fact6;

                count++;
             }
          }
       }
    }

    for(ii=0;ii<ROWS;ii++)
    {
      delete [] ders1[ii];
      delete [] ders2[ii];
    }
    delete [] ders1;
    delete [] ders2;
*/
    return;
}
