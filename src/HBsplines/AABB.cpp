

#include "AABB.h"


namespace myGeom
{

bool AABB::doIntersect(AABB&  bb)
{

#ifdef _DOMAIN2D
  if( (bb.minBB[0] > maxBB[0] ) || ( bb.maxBB[0] < minBB[0] ) ||
      (bb.minBB[1] > maxBB[1] ) || ( bb.maxBB[1] < minBB[1] ) )
  {
    return  false;
  }
  //
#else
  if( ( bb.minBB[0] > maxBB[0] ) || ( bb.maxBB[0] < minBB[0] ) ||
      ( bb.minBB[1] > maxBB[1] ) || ( bb.maxBB[1] < minBB[1] ) ||
      ( bb.minBB[2] > maxBB[2] ) || ( bb.maxBB[2] < minBB[2] )  )
  {
    return  false;
  }
#endif

  return  true;
}


bool AABB::within(myPoint&  pt)
{
#ifdef _DOMAIN2D
  if(pt[0] < minBB[0] || pt[0] > maxBB[0] ||
     pt[1] < minBB[1] || pt[1] > maxBB[1] )
  {
    return false;
  }
#else
  if(pt[0] < minBB[0] || pt[0] > maxBB[0] ||
     pt[1] < minBB[1] || pt[1] > maxBB[1] ||
     pt[2] < minBB[2] || pt[2] > maxBB[2])
  {
    return false;
  }
#endif

  return true;
}


}

