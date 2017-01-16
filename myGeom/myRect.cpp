

#include "myRect.h"

#include "BasisFunctionsLagrange.h"
#include "QuadratureUtil.h"

using namespace std;

namespace myGeom
{

myRect::~myRect()
{
}


double  myRect::volume()
{
    double val = ptList[0][0]*(ptList[1][1]-ptList[2][1]);
           val += ptList[1][0]*(ptList[2][1]-ptList[0][1]);
           val += ptList[2][0]*(ptList[0][1]-ptList[1][1]);
           val *= 0.5;
    return val;
}

int myRect::orientation()
{
  return POLY_ORIENT_RANDOM;
}



void myRect::reverseOrientation()
{ 
  return ; 
}



double  myRect::volumeFromGaussPoints(int nGP)
{
  vector<double>  gps1(nGP), gps2(nGP), gws(nGP);

  getGaussPointsTriangle(nGP, gps1, gps2, gws);
  
  int  ii, gp, nlbf=nVert, deg=1;
  double  Jac, detinv, xx, yy, dvol;
  double  N[4], dN_du[2][4];
  double  B[2][2], Binv[2][2];
  
  dvol = 0.0;
  for(gp=0; gp<nGP; gp++)
  {
      LagrangeBasisFunsTria(deg, gps1[gp], gps2[gp], N, dN_du[0], dN_du[1]);

      B[0][0] = B[0][1] = B[1][0] = B[1][1] = 0.0;
      for(ii=0; ii<nlbf; ii++)
      {
          xx = ptList[ii][0];
          yy = ptList[ii][1];

          B[0][0] +=  (xx * dN_du[0][ii]) ;
          B[1][0] +=  (xx * dN_du[1][ii]) ;
          B[0][1] +=  (yy * dN_du[0][ii]) ;
          B[1][1] +=  (yy * dN_du[1][ii]) ;
      }

      Jac = B[0][0]*B[1][1] - B[0][1]*B[1][0];
  
      //detinv = 1.0/Jac ;

      //Binv[0][0] =  B[1][1] * detinv;
      //Binv[0][1] = -B[0][1] * detinv;
      //Binv[1][0] = -B[1][0] * detinv;
      //Binv[1][1] =  B[0][0] * detinv;

      //for(ii=0; ii<nlbf; ii++)
      //{
        //dN_dx[ii] = dN_du[0][ii] * Binv[0][0] + dN_du[1][ii] * Binv[0][1];
        //dN_dy[ii] = dN_du[0][ii] * Binv[1][0] + dN_du[1][ii] * Binv[1][1];
      //}

      dvol += (gws[gp] * Jac);
  }

  return dvol;
}



void myRect::getGaussPointsCUTFEM(int nGP, vector<point3d>& gps, vector<double>& gws1 )
{
  vector<double>  gps1(nGP), gps2(nGP), gws(nGP);

  getGaussPointsTriangle(nGP, gps1, gps2, gws);
  
  int  ii, gp, nlbf=nVert, deg=1;
  double  Jac, detinv, xx, yy, x0, y0;
  double  N[4], dN_du[2][4];
  double  B[2][2], Binv[2][2];
  
  gps.resize(nGP);
  gws1.resize(nGP);
  
  for(gp=0; gp<nGP; gp++)
  {
      LagrangeBasisFunsTria(deg, gps1[gp], gps2[gp], N, dN_du[0], dN_du[1]);

      B[0][0] = B[0][1] = B[1][0] = B[1][1] = 0.0;
      x0 = y0 = 0.0;
      for(ii=0; ii<nlbf; ii++)
      {
          xx = ptList[ii][0];
          yy = ptList[ii][1];

          x0 += N[ii] * xx;
          y0 += N[ii] * yy;

          B[0][0] +=  (xx * dN_du[0][ii]) ;
          B[1][0] +=  (xx * dN_du[1][ii]) ;
          B[0][1] +=  (yy * dN_du[0][ii]) ;
          B[1][1] +=  (yy * dN_du[1][ii]) ;
      }

      Jac = B[0][0]*B[1][1] - B[0][1]*B[1][0];

      gps[gp][0] = x0;
      gps[gp][1] = y0;
      gws1[gp]   = gws[gp]*Jac;
  }

  return;
}


void myRect::computeNormal()
{
  // normal to a plane passing through 3 points defined by 
  // the conrners (P0, P1, P2) of a triangle
  // 
  // normal = (P1-P0) X (P2-P0)

  normal = (ptList[1]-ptList[0]).cross(ptList[2]-ptList[0]);

  normal.normalize();
}



double  myRect::distanceFromPoint(const myPoint& pt)
{
  // A plane with a point P0 and normal to the plane n
  // normal n has to be normalized
  // distance from a point P to that plane is
  // the project of (P-P0) on to the normal (n) to the plane
  //
  // distance = dotproduct(n, pt-p0)
  // 
  
  return  ( normal.dot(ptList[0] - pt) );  
}



int  myRect::IntersectWithRay(const Ray& ray1, vector<myPoint>& ptOut)
{
  //if( ((pt1[1]>ray1.Orig[1]) != (pt2[1]>ray1.Orig[1])) &&
    //(pt[0] < (ptList[j][0]-ptList[i][0]) * (pt[1]-ptList[i][1]) / (ptList[j][1]-ptList[i][1]) + ptList[i][0]) )

  // check if both the points of the line on the same side of the ray
  // If they lie on the same side then the ray does not intersect the line
  if( (ptList[0][1]>ray1.Orig[1]) == (ptList[1][1]>ray1.Orig[1]) )
  {
    return 0;
  }
  else
  {
    // If the points lie on the different sides of the ray then the ray intersects
    // the line segment
    return 1;
  }  
}



}












