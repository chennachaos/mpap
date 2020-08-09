

#ifndef incl_Intersection2D_h
#define incl_Intersection2D_h

#include "headersBasic.h"
#include "Ray.h"


inline int LineLineIntersection2D(myLine& ln1, myLine& ln2, myPoint& pout)
{
  // check if the bounding boxes of the two lines intersect
  
  double  r1, r2, l1, l2, t1, t2, b1, b2;
  
  r1 = max(ln1.ptList[0][0], ln1.ptList[1][0]);
  l1 = min(ln1.ptList[0][0], ln1.ptList[1][0]);
  t1 = max(ln1.ptList[0][1], ln1.ptList[1][1]);
  b1 = min(ln1.ptList[0][1], ln1.ptList[1][1]);

  r2 = max(ln2.ptList[0][0], ln2.ptList[1][0]);
  l2 = min(ln2.ptList[0][0], ln2.ptList[1][0]);
  t2 = max(ln2.ptList[0][1], ln2.ptList[1][1]);
  b2 = min(ln2.ptList[0][1], ln2.ptList[1][1]);

  //return !(l2 > r1 || r2 < l1 || t2 < b1 || b2 > t1);
  if( l2 > r1 )
  {
    return -1;
  }
  else
  {
    if( r2 < l1 )
      return -1;
    else
    {
      if( t2 < b1 )
        return -1;
      else
        if( b2 > t1 )
          return -1;
    }
  }

  myPoint  p1 = ln1.ptList[1] - ln1.ptList[0];
  myPoint  p2 = ln2.ptList[1] - ln2.ptList[0];
  
  double  m1 = p1[1]/p1[0];
  double  m2 = p2[1]/p2[0];
  
  double  tol = 1.0e-10;
  
  // check if the two lines are parallel by comparing the slopes
  // if they are parallel then they do not intersect
  if( (m2-m1) < tol )
    return -1;
  
  // if the lines are not parallel then they must intersect
  // find the intersection point
  
  double  c1 = ln1.ptList[0][1] - m1 * ln1.ptList[0][0];
  double  c2 = ln2.ptList[0][1] - m2 * ln2.ptList[0][0];
  pout[0] = (c1-c2)/(m2-m1);
  pout[1] = m1*pout[0] + c1;
  
  // check if the intersection point lies outside of any of the two line segments
  
  myPoint  p3 = pout - ln1.ptList[0];
 
  double param = p3.dot(p1)/p1.dot(p1);
  
  if (param < 0.0 || param > 1.0)
    return -1;

  p3 = pout - ln2.ptList[0];
  
  param = p3.dot(p2)/p2.dot(p2);

  if (param < 0.0 || param > 1.0)
    return -1;  
}



#endif