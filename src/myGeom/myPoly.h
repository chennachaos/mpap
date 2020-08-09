#ifndef INCL_myPoly_H
#define INCL_myPoly_H

#include "myGeomUtilities.h"
#include "headersBasic.h"
#include "AABB.h"
#include "Ray.h"

namespace myGeom
{

enum polyOrientation { POLY_ORIENT_CW=-1, POLY_ORIENT_CCW=1, POLY_ORIENT_RANDOM=0 };



class  myPoly
{
    //typedef  Vector3d  point_type;
    //typedef myPoint  point_type;

  protected:

    AABB  bbox;

    vector<myPoint>  ptList;
    myPoint  normal;
    int nVert, domainNum;

  public:

    myPoly();

    myPoly(vector<myPoint>&   pts1)
    {
      ptList = pts1;
      nVert = ptList.size();
    }

    virtual ~myPoly(){}

    void append(myPoint& pt1)
    {
      ptList.push_back( pt1 );
      nVert++;
    }

    myPoint&  GetPoint(int ii)
    {
      assert(ii < nVert);

      return  ptList[ii];
    }
    
    void  updatePoint(int ii, const myPoint& ptTemp)
    {
      assert(ii < nVert);

      ptList[ii] = ptTemp;
      return;
    }

    int GetNumVertices()
    {
      return nVert;
    }

    virtual  void  computeAABB();
    
    virtual AABB&  getAABB()
    { return bbox; }

    void  centroid(myPoint& cent);

    void SetDomainNumber(int dd)
    { domainNum = dd;      }

    int getDomainNumber()
    { return domainNum ;   }

    virtual void  computeNormal();
    
    myPoint&  GetNormal()
    { return  normal; }

    virtual double area();

    virtual double surfaceArea();

    virtual double volume();

    virtual double volumeFromGaussPoints(int deg);

    virtual void getGaussPointsCUTFEM(int deg, vector<myPoint>& gps, vector<double>& gws );

    virtual int orientation();

    virtual void reverseOrientation();

    virtual  double  distanceFromPoint(myPoint& pt);

    virtual  int  doIntersect(AABB& bb2);

    virtual  int  IntersectWithRay(Ray& ray1, vector<myPoint>& ptOut)
    {    return 0;    };

    virtual  int  IntersectWithLine(myLine& ln, vector<myPoint>& ptOut)
    {
      cout << " 'IntersectWithLine' is not defined .... " << endl;
      return 0;
    };

    virtual  int  IntersectWithTriangle(myTria& tria, vector<myPoint>& ptOut)
    {    return 0;    };

    virtual  int  IntersectWithRectangle(myRect& rect, vector<myPoint>& ptOut)
    {    return 0;    };

    virtual  bool  within(myPoint& pt);

    virtual  int  computeGeometry(myPoint&  param, myPoint&  _geom)
    { return -1;}

    virtual  int  computeNormal(myPoint&  param, myPoint&  _normal)
    { return -1;}

    virtual  int  computeBasisFunctions(myPoint& param, myPoint& geom, VectorXd&  Nb, double& Jac)
    { return -1;}
};

}


#endif //INCL_myPoly_H
