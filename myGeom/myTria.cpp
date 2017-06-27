

#include "myTria.h"

#include "BasisFunctionsLagrange.h"
#include "QuadratureUtil.h"

using namespace std;

namespace myGeom
{

myTria::~myTria()
{
}


void  myTria::computeAABB()
{
  //bbox.initialize();
  
  bbox.minBB[0] = min(min(ptList[0][0], ptList[1][0]), ptList[2][0]);
  bbox.maxBB[0] = max(max(ptList[0][0], ptList[1][0]), ptList[2][0]);

  bbox.minBB[1] = min(min(ptList[0][1], ptList[1][1]), ptList[2][1]);
  bbox.maxBB[1] = max(max(ptList[0][1], ptList[1][1]), ptList[2][1]);

  bbox.minBB[2] = min(min(ptList[0][2], ptList[1][2]), ptList[2][2]);
  bbox.maxBB[2] = max(max(ptList[0][2], ptList[1][2]), ptList[2][2]);

  yMin = bbox.minBB[1];
  yMax = bbox.maxBB[1];

  return;
}  


double  myTria::volume()
{
    double val = ptList[0][0]*(ptList[1][1]-ptList[2][1]);
           val += ptList[1][0]*(ptList[2][1]-ptList[0][1]);
           val += ptList[2][0]*(ptList[0][1]-ptList[1][1]);
           val *= 0.5;

    return val;
}

int myTria::orientation()
{
  return POLY_ORIENT_RANDOM;
}



void myTria::reverseOrientation()
{ 
  return ;
}



double  myTria::volumeFromGaussPoints(int nGP)
{
  // compute area of the quadrilateral using numerical integration
  // with Guass points
  
  vector<double>  gps1(nGP), gps2(nGP), gws(nGP);

  getGaussPointsTriangle(nGP, gps1, gps2, gws);

  int  ii, gp, nlbf=nVert, deg=1;
  double  Jac, detinv, xx, yy, dvol;
  double  N[3], dN_du[2][3];
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



void myTria::getGaussPointsCUTFEM(int nGP, vector<point3d>& gps, vector<double>& gws1 )
{
  vector<double>  gps1(nGP), gps2(nGP), gws(nGP);

  getGaussPointsTriangle(nGP, gps1, gps2, gws);
  
  int  ii, gp, nlbf=nVert, deg=1;
  double  Jac, detinv, xx, yy, x0, y0;
  double  N[3], dN_du[2][3];
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


double  myTria::distanceFromPoint(myPoint& pt)
{
  // A plane with a point P0 and normal to the plane n
  // normal n has to be normalized
  // distance from a point P to that plane is
  // the project of (P-P0) on to the normal (n) to the plane
  //
  // distance = dotproduct(n, pt-p0)
  // 

  return  ( normal.dot(pt) + d0 );
}



/*
int  myTria::IntersectWithRay(const Ray& ray1, vector<myPoint>& ptOut)
{
   // ray: O+tD ---- O-origin, D-direction
   // plane : ax+by+cz=d0

  // 1.) find the intersection of the ray with the plane of the triangle
  // 2.) Determine if that point lies inside or outside of the triangle
  
  // check if the ray intersects the plane of the triangle
  
  double  tmp = normal.dot(ray1.Dir);
  
  double  EPS = 1.0e-10;

  if( tmp > -EPS && tmp < EPS )
    return 0;

  double  t = -(d0 + ray1.Orig.dot(normal))/tmp;

  // if the intersection is to the negative side of the ray then
  // the intersection does not count

  if(t < EPS)
    return 0;
    
  // find the intersection point of the ray with the plane

  myPoint  ptInter = ray1.Orig + t*ray1.Dir;

  myPoint  E1 = ptList[1]-ptList[0];
  myPoint  E2 = ptList[2]-ptList[0];

  myPoint  ptTemp = ptInter-ptList[0];

  tmp = E1[0]*E2[1]-E1[1]*E2[0];

  double  u = (ptTemp[0]*E2[1]-ptTemp[1]*E2[0])/tmp;

  if(u < 0.0 || u > 1.0)
    return  0;

  double  v = (ptTemp[1]*E1[0]-ptTemp[0]*E1[1])/tmp;

  if(v < 0.0 || v > 1.0)
    return  0;

  if( (u+v) > 1.0 )
    return 0;

  return 1;
  //ptOut.push_back(ptTemp);
}
*/



//
int  myTria::IntersectWithRay(Ray& ray1, vector<myPoint>& ptOut)
{
   // ray: O+tD ---- O-origin, D-direction
   // plane : ax+by+cz=d0

  // 1.) find the intersection of the ray with the plane of the triangle
  // 2.) Determine if that point lies inside or outside of the triangle
  
  // check if the ray intersects the plane of the triangle

  myPoint  E1 = ptList[2]-ptList[0];
  myPoint  E2 = ptList[1]-ptList[0];
  myPoint   P = ray1.Dir.cross(E2);

  double  tmp = P.dot(E1);
  
  double  EPS = 1.0e-10;

  if( tmp > -EPS && tmp < EPS )
    return 0;

  tmp = 1.0/tmp;

  myPoint  T = ray1.Orig - ptList[0];
  myPoint  Q = T.cross(E1);


  double  u = tmp * T.dot(P);
  double  v = tmp * ray1.Dir.dot(Q);
  double  t = tmp * E2.dot(Q);

  //cout << " uvt = " << u << '\t' << v << '\t' << t << endl;

  if(u < 0.0 || u > 1.0)
    return 0;

  if(v < 0.0 || v > 1.0)
    return 0;
  
  if( (u+v) > 1.0 )
    return 0;

  // find the intersection point of the ray with the plane
  
  // if the intersection is to the negative side of the ray then
  // the intersection does not count

  if(t < EPS)
    return 0;

  return 1;

}
//





int  myTria::IntersectWithLine(myLine& ln, vector<myPoint>& ptOut)
{
  return 1;
}

/*
int  myTria::IntersectWithLine(const myPoint& P1, const myPoint& P2, vector<myPoint>& ptOut)
{
   // line: P1+t(P2-P1) ---- P1, P2 end points of the line
   // plane : ax+by+cz+d0=0

  // 1.) find the intersection of the line segment with the plane of the triangle
  // 2.) Determine if that point lies inside or outside of the triangle
  
  // check if the line segment intersects the plane of the triangle
  
  myPoint  Dir = P2-P1;
  
  double  tmp = normal.dot(Dir);
  
  double  EPS = 1.0e-7;

  if( tmp > -EPS && tmp < EPS )
    return 0;

  double  t = -(d0 + P1.dot(normal))/tmp;

  // if t is not in the range [0,1] then the intersection is 
  // outside the line segment and does not count

  if(t < EPS || t > 1.0)
    return 0;
    
  // find the intersection point of the ray with the plane

  myPoint  ptInter = P1 + t*Dir;

  myPoint  E1 = ptList[1]-ptList[0];
  myPoint  E2 = ptList[2]-ptList[0];

  myPoint  ptTemp = ptInter-ptList[0];

  tmp = E1[0]*E2[1]-E1[1]*E2[0];

  double  u = (ptTemp[0]*E2[1]-ptTemp[1]*E2[0])/tmp;

  if(u < 0.0 || u > 1.0)
    return  0;

  double  v = (ptTemp[1]*E1[0]-ptTemp[0]*E1[1])/tmp;

  if(v < 0.0 || v > 1.0)
    return  0;

  if( (u+v) > 1.0 )
    return 0;
  
  cout << " parameters         " << tmp << '\t' << t << '\t' << u  << '\t' << v << endl;
  cout << " intersection point " << ptInter[0] << '\t' << ptInter[1] << '\t' << ptInter[2] << endl;

  ptOut.push_back(ptInter);

  return 1;
}
*/


int  myTria::IntersectWithLine(myPoint& P1, const myPoint& P2, vector<myPoint>& ptOut)
{
   // line: P1+t(P2-P1) ---- P1, P2 end points of the line
   // plane : ax+by+cz+d0=0

  // 1.) find the intersection of the line segment with the plane of the triangle
  // 2.) Determine if that point lies inside or outside of the triangle
  
  // check if the line segment intersects the plane of the triangle
  
  myPoint  Dir = P2-P1;

  myPoint  E1 = ptList[2]-ptList[0];
  myPoint  E2 = ptList[1]-ptList[0];
  myPoint   P = Dir.cross(E2);

  double  tmp = P.dot(E1);
  
  double  EPS = 1.0e-7;

  if( tmp > -EPS && tmp < EPS )
    return 0;

  tmp = 1.0/tmp;

  myPoint  T = P1 - ptList[0];
  myPoint  Q = T.cross(E1);


  double  u = tmp * T.dot(P);
  double  v = tmp * Dir.dot(Q);
  double  t = tmp * E2.dot(Q);

  //cout << " uvt = " << u << '\t' << v << '\t' << t << endl;

  if(u < 0.0 || u > 1.0)
    return 0;

  if(v < 0.0 || v > 1.0)
    return 0;
  
  if( (u+v) > 1.0 )
    return 0;

  // find the intersection point of the ray with the plane
  
  // if t is not in the range [0,1] then the intersection is 
  // outside the line segment and does not count

  if(t < EPS || t > 1.0)
    return 0;

  myPoint  ptInter = P1 + t*Dir;

  //cout << " parameters         " << tmp << '\t' << t << '\t' << u  << '\t' << v << endl;
  //cout << " intersection point " << ptInter[0] << '\t' << ptInter[1] << '\t' << ptInter[2] << endl;

  ptOut.push_back(ptInter);

  return 1;
}
//


int  myTria::IntersectWithTriangle(myTria& triaIn, vector<myPoint>& ptOut)
{
  // triangle #1 with normal n1 and point P0
  // triangle #2 with normal n2 and point Q0

  // check if the bounding boxes of the two triangles intersect
  
  if( !bbox.doIntersect(triaIn.getAABB()) )
  {
    //cout << " bboxes do not intersect " << endl;
    return 0;
  }
  
  // if the bounding boxes intersect then perform further checks

  // check if the planes of the two triangles intersect
  // if the planes do not intersect then they are parallel
  //  // check if the two planes are collinear
  //  // if the planes are not collinear then there is no intersection
  //  // if the planes are collinear then check for intersection
  // if the planes intersect then the triangles may intersect
  // find the intersection line and check for containment
  
  myPoint  NN = normal.cross(triaIn.GetNormal());

  double  EPS = 1.0e-4;

  if( NN.squaredNorm() < EPS )
  {
    // the planes of both the triangles are parallel
    //cout << " parallel triangles " << endl;
    return 0;
  }
  else
  {
    // planes of both the triangles intersect
    
    // intersection of edges of triangle #1 with the triangle #2
    //cout << " hhhhhhhhhh " << endl;
    int  dd=0;
    //if( this->IntersectWithLine(tria.ptList[0], tria.ptList[1], ptOut) ||
      //  this->IntersectWithLine(tria.ptList[1], tria.ptList[2], ptOut) ||
      //  this->IntersectWithLine(tria.ptList[2], tria.ptList[0], ptOut) )
      //dd = 1;

    if( IntersectWithLine(triaIn.GetPoint(0), triaIn.GetPoint(1), ptOut) ||
        IntersectWithLine(triaIn.GetPoint(1), triaIn.GetPoint(2), ptOut) ||
        IntersectWithLine(triaIn.GetPoint(2), triaIn.GetPoint(0), ptOut) )
      dd = 1;

    // intersection of edges of triangle #2 with the triangle #1

    if( triaIn.IntersectWithLine(ptList[0], ptList[1], ptOut) ||
        triaIn.IntersectWithLine(ptList[1], ptList[2], ptOut) ||
        triaIn.IntersectWithLine(ptList[2], ptList[0], ptOut) )
      dd = 1;

    return dd;
  }
}


int  myTria::IntersectWithRectangle(myRect& rect, vector<myPoint>& ptOut)
{
  return 1;
}





int  myTria::computeGeometry(myPoint&  param, myPoint&  _geom)
{
  _geom = ptList[0];

  return 1;
}



void myTria::computeNormal()
{
  // normal to a plane passing through 3 points defined by 
  // the conrners (P0, P1, P2) of a triangle
  // 
  // normal = (P1-P0) X (P2-P0)

  normal = (ptList[2]-ptList[0]).cross(ptList[1]-ptList[0]);

  normal.normalize();
  
  d0 = -normal.dot(ptList[0]);
  
  return;
}




int  myTria::computeNormal(myPoint&  param, myPoint&  _normal)
{
  _normal = normal;

  return 1;
}




int  myTria::computeBasisFunctions(myPoint& param, myPoint& geom, VectorXd&  Nb, double& Jac)
{
    vector<double>  dN_dxi(nVert), dN_dzeta(nVert) ;

    LagrangeBasisFunsTria(1, param[0], param[1], &Nb[0], &dN_dxi[0], &dN_dzeta[0]);

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







