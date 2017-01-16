
#include "DAFunctionPlot.h"
#include "FunctionsEssGrp.h"
#include "Plot.h"
#include "MathBasic.h"


extern Plot plot;


using namespace std;




DAFunctionPlot::DAFunctionPlot(void)                       
{                                                  
  mxn = 20000;

  return;
}




	                                          
DAFunctionPlot::~DAFunctionPlot(void)                     
{         
  return;
}






void DAFunctionPlot::setup(int xx, int yy, int xx0, int yy0, 
                           int fc, int bc, int gc, int maxgl, bool kARatio, bool grid)
{
  DAWindow::setup(xx,yy,xx0,yy0,1.,1.,true);

  fCol = fc; 
  bCol = bc; 
  gCol = gc;

  mxgl = maxgl;

  kAR  = kARatio;

  showGrid = grid;

  return;
}








void DAFunctionPlot::draw(void)
{
  if (x.n < 1) return;

  int i, j, x1, y1, x2, y2;

  double dX  = mxx - mnx, 
         dY  = mxy - mny,
         mnX = mnx,
         mnY = mny,
         dgX, dgY, dd;

  if (dX < 1.e-10) { dX = 2.; mnX -= 1.; }
  if (dY < 1.e-10) { dY = 2.; mnY -= 1.; }

  // calculate grid spacing

  dd = dX * max(1.,double(dy)/double(dx)) * 10. / double(mxgl);
  i = 0; 
  while (dd <= 1.) { i++; dd *= 10.; } 
  while (dd > 10.) { i--; dd *= .1;  }
  if      (dd < 2.5) dgX = .2;
  else if (dd < 5.5) dgX = .5;
  else               dgX = 1.;
  for (j=0; j<i; j++) dgX *= .1;
  for (j=0; j>i; j--) dgX *= 10.;
  mnX = floor(mnX/dgX) * dgX;
  dX = dgX; while (dX+mnX < mxx-1.e-10) dX += dgX;

  dd = dY * max(1.,double(dx)/double(dy)) * 10. / double(mxgl); 
  i = 0; 
  while (dd <= 1.) { i++; dd *= 10.; } 
  while (dd > 10.) { i--; dd *= .1;  }
  if      (dd < 2.5) dgY = .2;
  else if (dd < 5.5) dgY = .5;
  else               dgY = 1.;
  for (j=0; j<i; j++) dgY *= .1;
  for (j=0; j>i; j--) dgY *= 10.;
  mnY = floor(mnY/dgY) * dgY;
  dY = dgY; while (dY+mnY < mxy-1.e-10) dY += dgY;

  // draw grid and put labels

  wipe(bCol);

  plot.setColour(gCol);

  dd = mnX;
  y1 = y0;
  y2 = y0+dy;
  while (dd < mnX+dX+dgX*.1)
  {
    x1 = x0 + roundToInt((dd-mnX)/dX*(double)dx);
    dd += dgX;
    essGrpDrawLine(x1,y1,x1,y2);
  }
  dd = mnY;
  x1 = x0;
  x2 = x0+dx;
  while (dd < mnY+dY+dgY*.1)
  {
    y1 = y0+dy - roundToInt((dd-mnY)/dY*(double)dy);
    dd += dgY;
    essGrpDrawLine(x1,y1,x2,y1);
  }
    
  plot.setColour(fCol);

  // draw curve

  x1 = x0+roundToInt((x[0]-mnX)/dX*(double)dx);
  y1 = y0+dy-roundToInt((y[0]-mnY)/dY*(double)dy);

  for (i=1; i<x.n; i++)
  {
    x2 = x0+roundToInt((x[i]-mnX)/dX*(double)dx);
    y2 = y0+dy-roundToInt((y[i]-mnY)/dY*(double)dy);

    essGrpDrawLine(x1,y1,x2,y2);
 
    x1 = x2;
    y1 = y2;
  }

  return;
}





void DAFunctionPlot::addXY(double xval, double yval)
{
  if (x.n == 0)
  {
    mnx = xval;
    mny = yval;

    mxx = xval;
    mxy = yval;
  }
  else
  {
    mnx = min(mnx,xval);
    mny = min(mny,yval);

    mxx = max(mxx,xval);
    mxy = max(mxy,yval);
  }
  x.append(xval);
  y.append(yval);

  if (x.n > mxn)
  {
    x.del(0); 
    y.del(0);
  }

  return;
}

