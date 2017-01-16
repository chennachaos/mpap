
#ifndef incl_myLine_h
#define incl_myLine_h

#include "headersBasic.h"
#include "Ray.h"
#include "myPoly.h"


namespace myGeom
{

class myLine : public myPoly
{
  protected:
    double  slope;
  
  public:
    
    myLine()
    {
      ptList.clear();
      ptList.push_back(myPoint(0.0,0.0,0.0));
      ptList.push_back(myPoint(1.0,1.0,1.0));
      nVert = 2;
    }

    myLine(myPoint& pt1, myPoint& pt2)
    {
      ptList.clear();
      ptList.push_back(pt1);
      ptList.push_back(pt2);
      nVert = 2;
    }

    myLine(vector<myPoint>&   pts1)
    {    
      assert(pts1.size() == 2);

      ptList = pts1;
      nVert = 2;
    }

    virtual ~myLine();
    
    virtual  void  computeNormal()
    {
      normal[0] =  ptList[1][1] - ptList[0][1];
      normal[1] = -ptList[1][0] + ptList[0][0];

      slope = -normal[0]/normal[1];

      normal.normalize();
    }
    
    void  computeSlope()
    {
      slope = (ptList[1][1] - ptList[0][1])/(ptList[1][0] - ptList[0][0]);
    }
    
    double  GetSlope()
    {  return  slope; }
    
    double  GetLength()
    {  return  ( (ptList[1]-ptList[0]).norm() ); }

    virtual  void  computeAABB();

    virtual  double  distanceFromPoint(myPoint& pt);

    virtual  int  IntersectWithRay(Ray& ray1, vector<myPoint>& ptOut);

    virtual  int  IntersectWithLine(myLine& ln, vector<myPoint>& ptOut);

    virtual  int  computeGeometry(myPoint&  param, myPoint&  _geom);

    virtual  int  computeNormal(myPoint&  param, myPoint&  _normal);

    virtual  int  computeBasisFunctions(myPoint& param, myPoint& geom, VectorXd&  Nb, double& Jac);

};

}

#endif

