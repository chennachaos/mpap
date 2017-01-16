
#ifndef INCL_MY_GEOM_UTILITIES_H
#define INCL_MY_GEOM_UTILITIES_H


//#include "headersBoost.h"
#include "headersBasic.h"
#include "util.h"
#include "../lib/DomainTree.h"

//#include <Eigen/Dense>

//namespace bg = boost::geometry;

using namespace std;
//using namespace Eigen;

//typedef Vector3d myPoint;

namespace  myGeom
{

class  myPoly;
class  myLine;
class  myTria;
class  myRect;
class  myTria;
class  myHex;

template<typename T>
bool pointExists(const vector<T>& ptList, T& pt1)
{
  double  EPS=1.0e-10;
  //T  ptTemp;
  bool val=false;
  for(int ii=0; ii<ptList.size(); ii++)
  {
    //ptTemp = ptList[ii] - pt1;
    //cout << ii << '\t' << ptTemp.squaredNorm() << endl;
    if( (ptList[ii] - pt1).squaredNorm() <= EPS )
    {
      val = true;
      break;
    }
  }

  return val;
}




inline double DistanceToLine2D(myPoint& pt1, myPoint& pt2, const myPoint& pt)
{
  myPoint p1 = pt1 - pt2;
  
  p1.normalize();
  
  return (p1(0)*(pt1(1)-pt(1)) - p1(1)*(pt1(0)-pt(0)));
}


inline int  myFindInt(vector<int>& data, int val)
{
  int index=-1;
  for(int ii=0; ii<data.size(); ii++)
  {
    if(data[ii] == val)
    {
      index = ii;
      break;
    }
  }
  return index;
}





/*
class  myPolygon
{
    typedef bg::model::d2::point_xy<double> point_type;

    typedef bg::model::polygon<point_type, false, false > polygontype; // closed ccw orientation

    private:

      //vector<myPoint>  pts;

      polygontype  poly;

    public:

      myPolygon(){}

      myPolygon(vector<myPoint>&   pts1)
      {
        for(int ii=0; ii<pts1.size(); ii++)
          poly.outer().push_back( point_type(pts1[ii][0], pts1[ii][1]) );
      }

      ~myPolygon(){}

      void append(myPoint>& pts1)
      {
          poly.outer().push_back( point_type(pts1[0], pts1[1]) );
      }

      void append(vector<myPoint>&   pts1)
      {
        for(int ii=0; ii<pts1.size(); ii++)
          poly.outer().push_back( point_type(pts1[ii][0], pts1[ii][1]) );
      }

      bool inside(myPoint& pt)
      {
        return (distance(pt) < 0.0);
      }

      virtual double distance(point_type& pt)
      {
        myPoint  pt1;
        pt1(0) = bg::get<0>(pt);
        pt1(1) = bg::get<1>(pt);

        return distance(pt1);
      }

      virtual double distance(myPoint& pt)
      {
        //vector<double>  dist(pts.size());
       
        VectorXd   dist(pts.size());
        int ii, nn=pts.size();
       
        for(int ii=0;ii<nn-1;ii++)
          dist[ii] = DistantToLine(pts[ii], pts[ii+1], pt);

        dist[nn-1] = DistantToLine(pts[nn-1], pts[0], pt);

        return dist.minCoeff();
      }
};
*/





}


#endif // INCL_MY_GEOM_UTILITIES_H

