#ifndef incl_AABB_h
#define incl_AABB_h

#include "headersBasic.h"
#include "myGeomUtilities.h"
#include "myConstants.h"

namespace myGeom
{

class  AABB
{
  public:

  myPoint  minBB, maxBB;
  
  AABB()
  {
    minBB[0] = minBB[1] = minBB[2] = 1e10;
    maxBB[0] = maxBB[1] = maxBB[2] = -1e10;
  }

  AABB(double x0, double x1)
  {
    minBB[0] = x0;
    maxBB[0] = x1;
  }

  AABB(double x0, double x1, double y0, double y1)
  {
    minBB[0] = x0;
    maxBB[0] = x1;

    minBB[1] = y0;
    maxBB[1] = y1;
  }
  
  AABB(double x0, double x1, double y0, double y1, double z0, double z1)
  {
    minBB[0] = x0;
    maxBB[0] = x1;

    minBB[1] = y0;
    maxBB[1] = y1;

    minBB[2] = z0;
    maxBB[2] = z1;
  }
  
  virtual ~AABB(){}
  
  void  updateData(int dir, double v0, double v1)
  {
    minBB[dir] = v0;
    maxBB[dir] = v1;
  }

  void  initialize()
  {
    minBB[0] = minBB[1] = minBB[2] = 1.0e10;
    maxBB[0] = maxBB[1] = maxBB[2] = -1.0e10;
  }
  
  void printSelf()
  {
    printf("AABB min = %12.6f \t %12.6f \t %12.6f \n", minBB[0], minBB[1], minBB[2]);
    printf("AABB max = %12.6f \t %12.6f \t %12.6f \n", maxBB[0], maxBB[1], maxBB[2]);
  }

  bool doIntersect(AABB&  bb);

  bool within(myPoint&  pt);

};

}


#endif


