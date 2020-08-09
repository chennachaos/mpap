#ifndef INCL_myHex_H
#define INCL_myHex_H

#include "myGeomUtilities.h"
#include "myPoly.h"
#include "headersBasic.h"


namespace myGeom
{

class  myHex: public myPoly
{
   private:

   public:

    myHex(){}

    virtual ~myHex(){}

    myHex(myPoint& pt1, myPoint& pt2, myPoint& pt3, myPoint& pt4)
    {
      ptList.clear();

      ptList.push_back(pt1);
      ptList.push_back(pt2);
      ptList.push_back(pt3);
      ptList.push_back(pt4);
      nVert = 4;
    }

    myHex(vector<myPoint>&   pts1)
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
};

}


#endif //INCL_myHex_H
