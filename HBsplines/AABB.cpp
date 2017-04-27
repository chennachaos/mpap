

#include "AABB.h"




namespace myGeom
{

bool AABB::doIntersect(AABB&  bb)
{
/*
  if( (bb.minBB[0] > maxBB[0] ) || ( bb.maxBB[0] < minBB[0] ) || 
      ( bb.minBB[1] > maxBB[1] ) || ( bb.maxBB[1] < minBB[1] ) )
  {
    return  false;
  }
  else
    return  true;
*/
//
  if( ( bb.minBB[0] > maxBB[0] ) || ( bb.maxBB[0] < minBB[0] ) || 
      ( bb.minBB[1] > maxBB[1] ) || ( bb.maxBB[1] < minBB[1] ) ||
      ( bb.minBB[2] > maxBB[2] ) || ( bb.maxBB[2] < minBB[2] )  )
  {
    return  false;
  }
  else
    return  true;
//
}


bool AABB::within(myPoint&  pt)
{
/*
  if(pt[0] < minBB[0] || pt[0] > maxBB[0] || 
    pt[1] < minBB[1] || pt[1] > maxBB[1] )
  {
    return false;
  }
*/
//
  if(pt[0] < minBB[0] || pt[0] > maxBB[0] || 
     pt[1] < minBB[1] || pt[1] > maxBB[1] || 
     pt[2] < minBB[2] || pt[2] > maxBB[2])
  {
    return false;
  }
//
  return true;
}

  
}



