
#include "headersBasic.h"
#include "Ray.h"
#include "myLine.h"

#include "BasisFunctionsLagrange.h"


namespace myGeom
{

myLine::~myLine()
{
}

void  myLine::computeAABB()
{
  //bbox.initialize();
  
  bbox.minBB[0] = min(ptList[0][0], ptList[1][0]);
  bbox.maxBB[0] = max(ptList[0][0], ptList[1][0]);

  bbox.minBB[1] = min(ptList[0][1], ptList[1][1]);
  bbox.maxBB[1] = max(ptList[0][1], ptList[1][1]);
  
  return;
}  



int  myLine::IntersectWithRay(Ray& ray1, vector<myPoint>& ptOut)
{
  // check if both the points of the line on the same side of the ray
  // If they lie on the same side then the ray does not intersect the line

  // If the points lie on the different sides of the ray then 
  // the ray intersects the line segment
  // Only intersections to the positive direction of the ray are counted as positive
  // Check if the line segment is to the left side of the ray's origin.
  // if so, then such an intersection is not counted

  if( ((ptList[0][1]>ray1.Orig[1]) != (ptList[1][1]>ray1.Orig[1])) &&
    (ray1.Orig[0] < (ptList[1][0]-ptList[0][0]) * (ray1.Orig[1]-ptList[0][1]) / (ptList[1][1]-ptList[0][1]) + ptList[0][0]) )
    return 1;
  else
    return 0;

  //if( ((ptList[0][1]>ray1.Orig[1]) != (ptList[1][1]>ray1.Orig[1])) &&
    //(ray1.Orig[0] <  (ray1.Orig[1]-ptList[0][1])/slope + ptList[0][0]) )
    //return 1;
  //else
    //return 0;
}



double myLine::distanceFromPoint(myPoint& pt)
{
  myPoint p1 = ptList[0] - ptList[1];
  
  p1.normalize();
  
  return ( p1[0]*(ptList[0][1]-pt[1]) - p1[1]*(ptList[0][0]-pt[0]) );
}




inline  double  dot2D(myPoint& pt1, myPoint& pt2)
{
  return  (pt1[0]*pt2[0] + pt1[1]*pt2[1]);
}



inline  double  perp2D(myPoint& pt1, myPoint& pt2)
{
  return  (pt1[0]*pt2[1] - pt1[1]*pt2[0]);
}


int myLine::IntersectWithLine(myLine& lnTemp, vector<myPoint>& ptOut)
{
   /* algorith is taken from
   ** http://geomalgorithms.com/a05-_intersect-1.html 
   */

  // check if the bounding boxes of the two lines intersect
  // if the bounding boxes do not intersect then the two line segments do not intersect

  if( ! bbox.doIntersect(lnTemp.GetAABB()) )
    return 0;

  // if the bounding boxes of the two line segments intersect then find the intersection point
  // and do further checking 

  myPoint  d1 = ptList[1] - ptList[0];
  myPoint  d2 = lnTemp.ptList[1] - lnTemp.ptList[0];
  myPoint  w = ptList[0] - lnTemp.ptList[0];

  double  D = perp2D(d1, d2);
  double  t1, t2, tol=1.0e-8, llim=-tol, ulim=1.0+tol;

  // check if the two lines are parallel by comparing the slopes
  // if they are parallel then they do not intersect

  if( abs(D) < tol )
  {
    // check if they are not collinear
    if( (perp2D(d1, w) != 0.0) || (perp2D(d2, w) != 0.0) )
      return 0;
    
    // parallel and collinear line segments
    // get the overlap, if any
  
    myPoint  w2 = ptList[1] - lnTemp.ptList[0];
  
    //r1 = d2[0]*lnTemp.ptList[0][0] + d2[1]*lnTemp.ptList[0][1] ;

    //t1 = (w[0]*lnTemp.ptList[0][0] + w[1]*lnTemp.ptList[0][1] )/r1;
    //t2 = (pt2[0]*lnTemp.ptList[0][0] + pt2[1]*lnTemp.ptList[0][1] )/r1;

    if( d2[0] != 0.0 )
    {
      t1 = w[0]/d2[0];
      t2 = w2[0]/d2[0];
    }
    else
    {
      t1 = w[1]/d2[1];
      t2 = w2[1]/d2[1];
    }

    if( t2 < llim || t1 > ulim) // disjoint line segments
      return 0;
    
    if( t1 < llim )
      ptOut.push_back(lnTemp.ptList[0]);

    //if( t1 >= 0.0 && t1 <= 1.0)
    if( t1 <= ulim)
      ptOut.push_back(ptList[0]);

    if( t2 > ulim )
      ptOut.push_back(lnTemp.ptList[1]);

    //if( t2 >= 0.0 && t2 <= 1.0)
    if( t2 >= llim )
      ptOut.push_back(ptList[1]);

    return 1;
  }

  // the segments are skew and may intersect in a point
  // check if the intersection point lies outside of any of the two line segments
  // get the intersect parameter for segment #1

  t1 = perp2D(d2, w) / D;

  if( t1 < llim || t1 > ulim )
    return 0;

  // get the intersect parameter for segment #2

  t2 = perp2D(d1, w) / D;

  if( t2 < 0.0 || t2 > 1.0 )
    return 0;

  // find the intersection point

  w = lnTemp.ptList[0] + t2*d2;

  //cout << pt1[0] << '\t' << pt1[1] << endl;
  ptOut.push_back(w);

  return 1;
}



int  myLine::computeGeometry(myPoint&  param, myPoint&  _geom)
{
  _geom = ptList[0];

  return 1;
}




int  myLine::computeNormal(myPoint&  param, myPoint&  _normal)
{
  _normal.setZero();
  //_normal = normal;

  _normal[0] =  ptList[1][1] - ptList[0][1];
  _normal[1] = -ptList[1][0] + ptList[0][0];

  //_normal.normalize();
  _normal /= _normal.norm();

  return 1;
}




int  myLine::computeBasisFunctions(myPoint& param, myPoint& geom, VectorXd&  Nb, double& Jac)
{
    vector<double>  dN1(nVert);

    Lagrange_BasisFuns1D(nVert-1, param[0], &Nb[0], &dN1[0]);

    //cout << xx[0] << '\t' << xx[1] << endl;

    double dx{0.0}, dy{0.0};
    geom.setZero();
    for(int ii=0; ii<nVert; ii++)
    {
      geom[0] +=  (ptList[ii][0] * Nb[ii]);
      geom[1] +=  (ptList[ii][1] * Nb[ii]);
      
      dx +=  (ptList[ii][0] * dN1[ii]);
      dy +=  (ptList[ii][1] * dN1[ii]);
    }
    
    Jac = sqrt(dx*dx+dy*dy);

    //double  du_dx = 1.0/Jac ;

    // Compute derivatives of basis functions w.r.t physical coordinates
    //for(ii=0; ii<nVert; ii++)
      //dN_dx[ii] = dN1[ii] * du_dx ;

  return 1;
}






}










