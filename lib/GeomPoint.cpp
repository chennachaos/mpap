
#include "GeomPoint.h"
#include "Plot.h"


extern Plot plot;



GeomPoint::GeomPoint(void)
{
  return;
}




GeomPoint::GeomPoint(double *X)
{
  x[0] = X[0];
  x[1] = X[1];

  return;
}




GeomPoint::~GeomPoint()
{
  return;
}





void GeomPoint::draw(int pointId)
{
  double d = (plot.dAct[0] + plot.dAct[1]) * .0055;

  plot.point(x,d,pointId);
 
  return;
}


