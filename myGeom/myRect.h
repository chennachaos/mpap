#ifndef INCL_myRect_H
#define INCL_myRect_H

#include "myGeomUtilities.h"
#include "myPoly.h"
#include "headersBasic.h"


namespace myGeom
{

class  myRect: public myPoly
{
   private:

   public:

    myRect(){}

    virtual ~myRect();

    myRect(myPoint& pt1, myPoint& pt2, myPoint& pt3)
    {
      ptList.clear();
      //ptList.erase();
      ptList.push_back(pt1);
      ptList.push_back(pt2);
      ptList.push_back(pt3);
      nVert = 3;
    }

    myRect(vector<myPoint>&   pts1)
    {    
      assert(pts1.size() == 3);

      ptList = pts1;
      nVert = 3;
    }
    
    virtual void  computeNormal();

    virtual double volume();

    virtual double volumeFromGaussPoints(int deg);

    virtual void getGaussPointsCUTFEM(int deg, vector<myPoint>& gps, vector<double>& gws );

    virtual int orientation();

    virtual void reverseOrientation();

    virtual  double  distanceFromPoint(const myPoint& pt);

    virtual  int  IntersectWithRay(const Ray& ray1, vector<myPoint>& ptOut);

};

}


#endif //INCL_myRect_H
