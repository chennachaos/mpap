

#include "myPoly.h"

#include "BasisFunctionsLagrange.h"
#include "QuadratureUtil.h"

using namespace std;

namespace myGeom
{

myPoly::myPoly()
{
  nVert = 0;
}


void  myPoly::computeNormal()
{
  return;
}


double myPoly::area()
{
/*
// area2D_Polygon(): compute the area of a 2D polygon
//  Input:  int n = the number of vertices in the polygon
//          Point* V = an array of n+1 vertex points with V[n]=V[0]
//  Return: the (float) area of the polygon
float
area2D_Polygon( int n, Point* V )
{
    float area = 0;
    int  i, j, k;   // indices

    if (n < 3) return 0;  // a degenerate polygon

    for (i=1, j=2, k=0; i<n; i++, j++, k++) {
        area += V[i].x * (V[j].y - V[k].y);
    }
    area += V[n].x * (V[1].y - V[n-1].y);  // wrap-around term
    return area / 2.0;
}
*/
  return 0.0;
}


double myPoly::surfaceArea()
{
  return 0.0;
}


double myPoly::volume()
{
  return 0.0;
}


double myPoly::volumeFromGaussPoints(int deg)
{
  return 0.0;
}


void myPoly::getGaussPointsCUTFEM(int deg, vector<myPoint>& gps, vector<double>& gws )
{
  return;
}


int myPoly::orientation()
{
/*
// orientation2D_Polygon(): test the orientation of a simple 2D polygon
//  Input:  int n = the number of vertices in the polygon
//          Point* V = an array of n+1 vertex points with V[n]=V[0]
//  Return: >0 for counterclockwise 
//          =0 for none (degenerate)
//          <0 for clockwise
//  Note: this algorithm is faster than computing the signed area.
int
orientation2D_Polygon( int n, Point* V )
{
    // first find rightmost lowest vertex of the polygon
    int rmin = 0;
    int xmin = V[0].x;
    int ymin = V[0].y;

    for (int i=1; i<n; i++) {
        if (V[i].y > ymin)
            continue;
        if (V[i].y == ymin) {   // just as low
            if (V[i].x < xmin)  // and to left
                continue;
        }
        rmin = i;      // a new rightmost lowest vertex
        xmin = V[i].x;
        ymin = V[i].y;
    }

    // test orientation at the rmin vertex
    // ccw <=> the edge leaving V[rmin] is left of the entering edge
    if (rmin == 0)
        return isLeft( V[n-1], V[0], V[1] );
    else
        return isLeft( V[rmin-1], V[rmin], V[rmin+1] );
}
*/
  return 0;
}


void myPoly::reverseOrientation()
{
  return ;
}

//int  myPoly::IntersectWithRay(Ray& ray1, vector<myPoint>& ptOut)
//{  return 0;}


//int  myPoly::IntersectWithLine(myLine& ln, vector<myPoint>& ptOut)
//{  return 0;}


//int  myPoly::IntersectWithTriangle(myTria& tria, vector<myPoint>& ptOut)
//{  return 0;}


//int  myPoly::IntersectWithRectangle(myRect& rect, vector<myPoint>& ptOut)
//{  return 0;}


void  myPoly::computeAABB()
{
  //bbox.initialize();
  
  int ii, jj;

  bbox.minBB[0] = ptList[0][0] ;
  bbox.minBB[1] = ptList[0][1] ;
  bbox.minBB[2] = ptList[0][2] ;

  bbox.maxBB[0] = ptList[0][0] ;
  bbox.maxBB[1] = ptList[0][1] ;
  bbox.maxBB[2] = ptList[0][2] ;

  for(ii=1; ii<nVert; ii++)
  {
    for(jj=0; jj<3; jj++)
    {
      bbox.minBB[jj] = min(bbox.minBB[jj], ptList[ii][jj]);
      bbox.maxBB[jj] = max(bbox.maxBB[jj], ptList[ii][jj]);

      //if( ptList[ii][jj] < bbox.minBB[jj]  )
        //bbox.minBB[jj] = ptList[ii][jj] ;

      //if( ptList[ii][jj] > bbox.maxBB[jj]  )
        //bbox.maxBB[jj] = ptList[ii][jj] ;
    }
  }

  return;
}
  
  

bool myPoly::within(myPoint& pt)
{
  //return (distance(pt) < 0.0);
  
  // check if the point is outside/inside the boundingbox of the polygon
  // if the point is outside the boundingbox then
  // the point is outside the polygon
  // otherwise, the point may be inside the polygon and needs further checking

  if( bbox.within(pt) )
  {
    int i, j, c=0;
    for(i=0; i<(nVert-1); i++)
    {
      j = i+1;
      if ( ((ptList[i][1]>pt[1]) != (ptList[j][1]>pt[1])) &&
       (pt[0] < (ptList[j][0]-ptList[i][0]) * (pt[1]-ptList[i][1]) / (ptList[j][1]-ptList[i][1]) + ptList[i][0]) )
        c = !c;
    }

    return c;
  }
  else
  {
    return false;
  }  
  
}



void myPoly::centroid(myPoint& cent)
{
  cent.setZero();

  for(int ii=0;ii<nVert;ii++)
    cent += ptList[ii];

  cent /= nVert;

  return ;
}


double myPoly::distanceFromPoint(myPoint& pt)
{
  VectorXd   dist(nVert);

  for(int ii=0;ii<nVert-1;ii++)
    dist[ii] = DistanceToLine2D(ptList[ii], ptList[ii+1], pt);

  dist[nVert-1] = DistanceToLine2D(ptList[nVert-1], ptList[0], pt);

  return dist.minCoeff();
}


int  myPoly::doIntersect(AABB& bb2)
{
  return bbox.doIntersect(bb2);
}







}






