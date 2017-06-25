

#include "myQuad.h"

#include "BasisFunctionsLagrange.h"
#include "QuadratureUtil.h"

using namespace std;

namespace myGeom
{

myQuad::~myQuad()
{
}


double  myQuad::volume()
{
    double val = ptList[0][0]*(ptList[1][1]-ptList[2][1]);
           val += ptList[1][0]*(ptList[2][1]-ptList[0][1]);
           val += ptList[2][0]*(ptList[0][1]-ptList[1][1]);
           val *= 0.5;

    return val;
}



int myQuad::orientation()
{
  return POLY_ORIENT_RANDOM;
}



void myQuad::reverseOrientation()
{ 
  return ; 
}



double  myQuad::volumeFromGaussPoints(int nGP)
{
  // compute area of the quadrilateral using numerical integration
  // with Guass points

  vector<double>  gps1(nGP), gps2(nGP), gws(nGP);

  getGaussPointsQuad(nGP, gps1, gps2, gws);
  
  int  ii, gp, nlbf=nVert, deg=1;
  double  Jac, detinv, xx, yy, dvol;
  double  N[4], dN_du[2][4];
  double  B[2][2], Binv[2][2];

  dvol = 0.0;
  for(gp=0; gp<nGP; gp++)
  {
      LagrangeBasisFunsQuad(deg, gps1[gp], gps2[gp], N, dN_du[0], dN_du[1]);

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

      dvol += (gws[gp] * Jac);
  }

  return dvol;
}



void myQuad::getGaussPointsCUTFEM(int nGP, vector<point3d>& gps, vector<double>& gws1 )
{
  cout << " myQuad::getGaussPointsCUTFEM ... NOT implemented " << endl;

  return;
}



double  myQuad::distanceFromPoint(const myPoint& pt)
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



int  myQuad::IntersectWithRay(const Ray& ray1, vector<myPoint>& ptOut)
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


int  myQuad::computeGeometry(myPoint&  param, myPoint&  _geom)
{
  _geom = ptList[0];

  return 1;
}



void myQuad::computeNormal()
{
  // P2------P3
  // |       |
  // |       |
  // |       |
  // |       |
  // P0------P1
  //
  // normal to a plane passing through 4 points defined by 
  // the conrners (P0, P1, P3, P2) of a quadrilateral
  // 
  // normal = (P1-P0) X (P2-P0)

  normal = (ptList[1]-ptList[0]).cross(ptList[2]-ptList[0]);

  normal.normalize();
}



int  myQuad::computeNormal(myPoint&  param, myPoint&  _normal)
{
  _normal = normal;

  return 1;
}




int  myQuad::computeBasisFunctions(myPoint& param, myPoint& geom, VectorXd&  Nb, double& Jac)
{
    vector<double>  dN_dxi(nVert), dN_dzeta(nVert) ;

    LagrangeBasisFunsQuad(1, param[0], param[1], &Nb[0], &dN_dxi[0], &dN_dzeta[0]);

    Vector3d  dx1, dx2;
    geom.setZero();
    dx1.setZero();
    dx2.setZero();

    for(int ii=0; ii<nVert; ii++)
    {
      geom[0] +=  (ptList[ii][0] * Nb[ii]);
      geom[1] +=  (ptList[ii][1] * Nb[ii]);
      geom[2] +=  (ptList[ii][2] * Nb[ii]);
      
      dx1(0) +=  (ptList[ii][0] * dN_dxi[ii]);
      dx2(0) +=  (ptList[ii][0] * dN_dzeta[ii]);

      dx1(1) +=  (ptList[ii][1] * dN_dxi[ii]);
      dx2(1) +=  (ptList[ii][1] * dN_dzeta[ii]);

      dx1(2) +=  (ptList[ii][2] * dN_dxi[ii]);
      dx2(2) +=  (ptList[ii][2] * dN_dzeta[ii]);
    }
    dx1 = dx1.cross(dx2);
    Jac = dx1.norm();

    //double  du_dx = 1.0/Jac ;

    // Compute derivatives of basis functions w.r.t physical coordinates
    //for(ii=0; ii<nVert; ii++)
      //dN_dx[ii] = dN1[ii] * du_dx ;

  return 1;
}



}










