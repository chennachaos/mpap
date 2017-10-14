#ifndef INCL_MYTET_H
#define INCL_MYTET_H

#include "myGeomUtilities.h"
#include "myPoly.h"
#include "headersBasic.h"


namespace myGeom
{

class  myTet: public myPoly
{
   private:

   public:

    myTet(){}

    virtual ~myTet();

    myTet(myPoint& pt1, myPoint& pt2, myPoint& pt3, myPoint& pt4)
    {
      ptList.clear();

      ptList.push_back(pt1);
      ptList.push_back(pt2);
      ptList.push_back(pt3);
      ptList.push_back(pt4);
      nVert = 4;
    }

    myTet(vector<myPoint>&   pts1)
    {    
      assert(pts1.size() == 4);

      ptList = pts1;
      nVert = 4;
    }

    virtual double surfaceArea();

    virtual double volume();

    virtual double volumeFromGaussPoints(int deg);

    virtual void getGaussPointsCUTFEM(int deg, vector<myPoint>& gps, vector<double>& gws );

    virtual int orientation();

    virtual void reverseOrientation();

    virtual  double  distanceFromPoint(const myPoint& pt);

    virtual  int  IntersectWithRay(const Ray& ray1, vector<myPoint>& ptOut);
};

}


#endif //INCL_myTet_H
