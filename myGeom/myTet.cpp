

#include "myTet.h"

#include "BasisFunctionsLagrange.h"
#include "QuadratureUtil.h"

using namespace std;

namespace myGeom
{

myTet::~myTet()
{
}

double  myTet::volume()
{
    MatrixXd  AA(4,4);
    
    for(int ii=0;ii<4; ii++)
    {
      for(int jj=0; jj<3; jj++)
        AA(ii, jj) = ptList[ii][jj];

      AA(ii, 3) = 1.0;
    }
    
    return (abs(AA.determinant()) / 6.0);
}

double myTet::surfaceArea()
{
  return 0.0;
}


int myTet::orientation()
{
  return POLY_ORIENT_RANDOM;
}



void myTet::reverseOrientation()
{ 
  return ; 
}



double  myTet::volumeFromGaussPoints(int nGP)
{
  vector<double>  gps1(nGP), gps2(nGP), gps3(nGP), gws(nGP);

  getGaussPointsTet(nGP, gps1, gps2, gps3, gws);
  
  int  ii, gp, nlbf=nVert, deg=1;
  double  xx, yy, zz, dvol;
  double  N[4], dN_du[3][4];
  MatrixXd  Bmat(3,3);
  
  dvol = 0.0;
  for(gp=0; gp<nGP; gp++)
  {
      LagrangeBasisFunsTet(deg, gps1[gp], gps2[gp], gps3[gp], N, dN_du[0], dN_du[1], dN_du[2]);

      Bmat.setZero();
      for(ii=0; ii<nlbf; ii++)
      {
          xx = ptList[ii][0];
          yy = ptList[ii][1];
          zz = ptList[ii][2];

          Bmat(0,0) +=  (xx * dN_du[0][ii]) ;
          Bmat(1,0) +=  (xx * dN_du[1][ii]) ;
          Bmat(2,0) +=  (xx * dN_du[2][ii]) ;
          
          Bmat(0,1) +=  (yy * dN_du[0][ii]) ;
          Bmat(1,1) +=  (yy * dN_du[1][ii]) ;
          Bmat(2,1) +=  (yy * dN_du[2][ii]) ;

          Bmat(0,2) +=  (zz * dN_du[0][ii]) ;
          Bmat(1,2) +=  (zz * dN_du[1][ii]) ;
          Bmat(2,2) +=  (zz * dN_du[2][ii]) ;
      }

      dvol += (gws[gp] * Bmat.determinant());
  }

  return dvol;
}



void myTet::getGaussPointsCUTFEM(int nGP, vector<point3d>& gps, vector<double>& gws1 )
{
  vector<double>  gps1(nGP), gps2(nGP), gps3(nGP), gws(nGP);

  getGaussPointsTet(nGP, gps1, gps2, gps3, gws);
  
  int  ii, gp, nlbf=nVert, deg=1;
  double  Jac, detinv, xx, yy, zz, x0, y0, z0;
  double  N[4], dN_du[3][4];
  MatrixXd  Bmat(3,3);
  
  gps.resize(nGP);
  gws1.resize(nGP);
  
  for(gp=0; gp<nGP; gp++)
  {
      LagrangeBasisFunsTet(deg, gps1[gp], gps2[gp], gps3[gp], N, dN_du[0], dN_du[1], dN_du[2]);

      Bmat.setZero();
      x0 = y0 = z0 = 0.0;
      for(ii=0; ii<nlbf; ii++)
      {
          xx = ptList[ii][0];
          yy = ptList[ii][1];
          zz = ptList[ii][2];

          x0 += N[ii] * xx;
          y0 += N[ii] * yy;
          z0 += N[ii] * zz;

          Bmat(0,0) +=  (xx * dN_du[0][ii]) ;
          Bmat(1,0) +=  (xx * dN_du[1][ii]) ;
          Bmat(2,0) +=  (xx * dN_du[2][ii]) ;
          
          Bmat(0,1) +=  (yy * dN_du[0][ii]) ;
          Bmat(1,1) +=  (yy * dN_du[1][ii]) ;
          Bmat(2,1) +=  (yy * dN_du[2][ii]) ;

          Bmat(0,2) +=  (zz * dN_du[0][ii]) ;
          Bmat(1,2) +=  (zz * dN_du[1][ii]) ;
          Bmat(2,2) +=  (zz * dN_du[2][ii]) ;
      }

      gps[gp][0] = x0;
      gps[gp][1] = y0;
      gps[gp][2] = z0;
      gws1[gp]   = gws[gp]*Bmat.determinant();
  }

  return;
}


double  myTet::distanceFromPoint(const myPoint& pt)
{
  return -11110.0;  
}



int  myTet::IntersectWithRay(const Ray& ray1, vector<myPoint>& ptOut)
{
  return 0;  
}
 
}


