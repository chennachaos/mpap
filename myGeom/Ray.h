
#ifndef incl_Ray_h
#define incl_Ray_h

//#include "headersBasic.h"
#include "util.h"

namespace  myGeom
{

class Ray
{
    public:
    
      myPoint Orig, Dir, normal;
      double  slope;
    
    Ray()
    {
      Orig.setZero();
      Dir = Orig;
      Dir[0] = 1.0;
    }
    
    Ray(myPoint& p1, myPoint& p2)
    {
      Orig = p1;
      Dir  = p2;
    }
    
    virtual ~Ray(){}
    
    void  updateOrigin(const myPoint& pt)
    {
      Orig = pt;
      return;
    }
    
    void  updateDirection(const myPoint& pt)
    {
      Dir = pt;
      computeNormal();
      computeSlope();

      return;
    }
    
    myPoint&  GetOrigin()
    {  return  Orig; }
    
    myPoint&  GetDirection()
    {  return  Dir; }
    
    void  computeNormal()
    {
      normal[0] =  Dir[1];
      normal[1] = -Dir[0];
    }

    myPoint&  GetNormal()
    { return  normal; }

    void  computeSlope()
    {
      slope = Dir[1]/Dir[0];
    }

    double  GetSlope()
    {  return  slope; }
    
    //myPoint

};

}


#endif

