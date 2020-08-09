#ifndef INCL_MYTRIA_H
#define INCL_MYTRIA_H

#include "myGeomUtilities.h"
#include "myPoly.h"
#include "headersBasic.h"


namespace myGeom
{

class  myTria: public myPoly
{
   private:
    double  d0, yMin, yMax;

   public:

    myTria(){}

    virtual ~myTria();

    myTria(myPoint& pt1, myPoint& pt2, myPoint& pt3)
    {
      ptList.clear();
      //ptList.erase();
      ptList.push_back(pt1);
      ptList.push_back(pt2);
      ptList.push_back(pt3);
      nVert = 3;
    }

    myTria(vector<myPoint>&   pts1)
    {    
      assert(pts1.size() == 3);

      ptList = pts1;
      nVert = 3;
    }
    
    virtual void  computeNormal();

    virtual  void  computeAABB();

    virtual double volume();

    virtual double volumeFromGaussPoints(int deg);

    virtual void getGaussPointsCUTFEM(int deg, vector<myPoint>& gps, vector<double>& gws );

    virtual int orientation();

    virtual void reverseOrientation();

    virtual  double  distanceFromPoint(myPoint& pt);

    virtual  int  IntersectWithRay(Ray& ray1, vector<myPoint>& ptOut);

    int  IntersectWithLine(myPoint& P1, const myPoint& P2, vector<myPoint>& ptOut);
    
    virtual  int  IntersectWithLine(myLine& ln, vector<myPoint>& ptOut);

    virtual  int  IntersectWithTriangle(myTria& tria, vector<myPoint>& ptOut);

    virtual  int  IntersectWithRectangle(myRect& rect, vector<myPoint>& ptOut);

    virtual  int  computeGeometry(myPoint&  param, myPoint&  _geom);

    virtual  int  computeNormal(myPoint&  param, myPoint&  _normal);

    virtual  int  computeBasisFunctions(myPoint& param, myPoint& geom, VectorXd&  Nb, double& Jac);

};

}


#endif //INCL_MYTRIA_H
