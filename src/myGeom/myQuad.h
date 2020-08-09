#ifndef INCL_myQuad_H
#define INCL_myQuad_H

#include "myGeomUtilities.h"
#include "myPoly.h"
#include "headersBasic.h"

/******************************
3--------2
|        |
|        |
|        |
0--------1
******************************/


namespace myGeom
{

class  myQuad: public myPoly
{
  private:

  public:

    myQuad(){}

    virtual ~myQuad();

    myQuad(myPoint& pt1, myPoint& pt2, myPoint& pt3, myPoint& pt4)
    {
      ptList.clear();
      //ptList.erase();
      ptList.push_back(pt1);
      ptList.push_back(pt2);
      ptList.push_back(pt3);
      ptList.push_back(pt4);
      nVert = 4;
    }

    myQuad(vector<myPoint>&   pts1)
    {    
      assert(pts1.size() == 4);

      ptList = pts1;
      nVert = 4;
    }
    
    virtual void  computeNormal();

    virtual double volume();

    virtual double volumeFromGaussPoints(int deg);

    virtual void getGaussPointsCUTFEM(int deg, vector<myPoint>& gps, vector<double>& gws );

    virtual int orientation();

    virtual void reverseOrientation();

    virtual  double  distanceFromPoint(const myPoint& pt);

    virtual  int  IntersectWithRay(const Ray& ray1, vector<myPoint>& ptOut);

    virtual  int  computeGeometry(myPoint&  param, myPoint&  _geom);

    virtual  int  computeNormal(myPoint&  param, myPoint&  _normal);

    virtual  int  computeBasisFunctions(myPoint& param, myPoint& geom, VectorXd&  Nb, double& Jac);

};

}


#endif //INCL_myQuad_H
