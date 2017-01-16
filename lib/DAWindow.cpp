
#include "DAWindow.h"
#include "Plot.h"


extern Plot plot;


using namespace std;




DAWindow::DAWindow(void)                       
{                                                  
  return;
}




	                                          
DAWindow::~DAWindow(void)                     
{         
  return;
}






void DAWindow::setup(int xx, int yy, int dxx, int dyy, double maxX, double maxY, bool kAR)
{
  x0 = xx;
  y0 = yy;
  dx = dxx;
  dy = dyy;

  mx = maxX;
  my = maxY;

  if (kAR) ; //....

  return;
}





void DAWindow::wipe(int c)
{
  plot.setColour(c);

  plot.iPolyPnt[0] = x0;
  plot.iPolyPnt[1] = y0;
  plot.iPolyPnt[2] = x0+dx;
  plot.iPolyPnt[3] = y0;
  plot.iPolyPnt[4] = x0+dx;
  plot.iPolyPnt[5] = y0+dy;
  plot.iPolyPnt[6] = x0;
  plot.iPolyPnt[7] = y0+dy;

  essGrpFillPoly(plot.iPolyPnt,4);
 
  return;
}




void DAWindow::frame(int c, int d)
{
  plot.setColour(c);

  plot.iPolyPnt[ 0] = x0;
  plot.iPolyPnt[ 1] = y0;
  plot.iPolyPnt[ 2] = x0+dx;
  plot.iPolyPnt[ 3] = y0;
  plot.iPolyPnt[ 4] = x0+dx;
  plot.iPolyPnt[ 5] = y0+dy;
  plot.iPolyPnt[ 6] = x0;
  plot.iPolyPnt[ 7] = y0+dy;
  plot.iPolyPnt[ 8] = x0;
  plot.iPolyPnt[ 9] = y0;
  plot.iPolyPnt[10] = x0+d;
  plot.iPolyPnt[11] = y0+d;
  plot.iPolyPnt[12] = x0+dx-d;
  plot.iPolyPnt[13] = y0+d;
  plot.iPolyPnt[14] = x0+dx-d;
  plot.iPolyPnt[15] = y0+dy-d;
  plot.iPolyPnt[16] = x0+d;
  plot.iPolyPnt[17] = y0+dy-d;
  plot.iPolyPnt[18] = x0+d;
  plot.iPolyPnt[19] = y0+d;

  essGrpFillPoly(plot.iPolyPnt,10);

  return;
}




void DAWindow::fillRectangle(int c, double xx, double yy, double dxx, double dyy, bool flag)
{
  plot.setColour(c); 

  int *PNT = plot.iPolyPnt;

  PNT[0] = x0 + roundToInt(      xx/mx * dx);
  PNT[1] = y0 + roundToInt(      yy/mx * dy);
  PNT[2] = x0 + roundToInt((xx+dxx)/mx * dx);
  PNT[3] = PNT[1];
  PNT[4] = PNT[2];
  PNT[5] = y0 + roundToInt((yy+dyy)/mx * dy);
  PNT[6] = PNT[0];
  PNT[7] = PNT[5];

  if (PNT[0] > PNT[2]) { PNT[0] = PNT[2]; PNT[2] = PNT[6]; PNT[4] = PNT[2]; PNT[6] = PNT[0]; }

  if (PNT[1] > PNT[5]) { PNT[1] = PNT[5]; PNT[5] = PNT[3]; PNT[3] = PNT[1]; PNT[7] = PNT[5]; }


  if (flag) // if rectangle is tiny, plot at least one pixel
  {
    if (PNT[0] == PNT[2]) { PNT[2]++; PNT[4]++; }
    if (PNT[3] == PNT[5]) { PNT[5]++; PNT[7]++; }
  }

  if (PNT[2] - PNT[0] > 8) { PNT[2]--; PNT[4]--; }
  if (PNT[5] - PNT[3] > 8) { PNT[5]--; PNT[7]--; }

  essGrpFillPoly(PNT,4);

  return;
}


