

#include "NodeGeomLink.h"






NodeGeomLink::NodeGeomLink(void)
{
  id = -1;

  return;
}








NodeGeomLink::NodeGeomLink(NodeGeomLink &ndGmLnk)
{
  for (int i=0; i<ndGmLnk.geomObj.n; i++) geomObj.append(ndGmLnk.geomObj[i]);

  id = ndGmLnk.id;

  return;
}








void NodeGeomLink::isOn(void *g0, void *g1, void *g2, void *g3)
{
  geomObj.append(g0);

  if (g1 != NULL) 
  {
    geomObj.append(g1);
    if (g2 != NULL) 
    {
      geomObj.append(g2);
      if (g3 != NULL) geomObj.append(g3);
    }
  }
  return;
}






void NodeGeomLink::isNotOn(void *g0, void *g1, void *g2, void *g3)
{
  int j, n = geomObj.n;

  j = 0; while (j < n && geomObj[j] != g0) j++;
  if (j < n) { geomObj.del(j); n--; }

  if (g1 != NULL) 
  {
    j = 0; while (j < n && geomObj[j] != g1) j++;
    if (j < n) { geomObj.del(j); n--; }
    
    if (g2 != NULL)
    {
      j = 0; while (j < n && geomObj[j] != g2) j++;
      if (j < n) { geomObj.del(j); n--; }

      if (g3 != NULL)
      {
        j = 0; while (j < n && geomObj[j] != g3) j++;
        if (j < n) geomObj.del(j); 
      }
    }
  }
  return;
}








bool NodeGeomLink::onThis(GeomObject *gO)
{
  int j = 0, n = geomObj.n;
  while (j < n && geomObj[j] != (void*)gO) j++;
  if (j < n) return true;
  return false;
}




bool NodeGeomLink::onThis(void *gO)
{
  int j = 0, n = geomObj.n;
  while (j < n && geomObj[j] != gO) j++;
  if (j < n) return true;
  return false;
}





bool NodeGeomLink::onSurface(void)
{
  int j = 0, n = geomObj.n;
  while (j < n && !((GeomObject*)(geomObj[j]))->isType(SURFACE)) j++;
  if (j < n) return true;
  return false;
}




GeomSurface* NodeGeomLink::whichSurface(void)
{
  int j = 0, n = geomObj.n;
  while (j < n && !((GeomObject*)(geomObj[j]))->isType(SURFACE)) j++;
  if (j < n) return (GeomSurface*)(geomObj[j]);
  return NULL;
}







bool NodeGeomLink::onSpline(void)
{
  int j = 0, n = geomObj.n;
  while (j < n && !((GeomObject*)(geomObj[j]))->isType(SPLINE)) j++;
  if (j < n) return true;
  return false;
}




GeomSpline* NodeGeomLink::whichSpline(void)
{
  int j = 0, n = geomObj.n;
  while (j < n && !((GeomObject*)(geomObj[j]))->isType(SPLINE)) j++;
  if (j < n) return (GeomSpline*)(geomObj[j]);
  return NULL;
}






bool NodeGeomLink::onPoint(void)
{
  int j = 0, n = geomObj.n;
  while (j < n && !((GeomObject*)(geomObj[j]))->isType(POINT)) j++;
  if (j < n) return true;
  return false;
}





GeomPoint* NodeGeomLink::whichPoint(void)
{
  int j = 0, n = geomObj.n;
  while (j < n && !((GeomObject*)(geomObj[j]))->isType(POINT)) j++;
  if (j < n) return (GeomPoint*)(geomObj[j]);
  return NULL;
}




