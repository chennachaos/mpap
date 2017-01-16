
#include <Xm/Xm.h>
#include <iostream>


#include "MathBasic.h"
#include "UnixGUI.h"


extern UnixGUI unixGUI;
extern bool noGUI;


void essGrpFillCircle(int *xx, int *d)
{
  if (noGUI) return;
  
  int a1 = 0, a2 = 23040, 
      x = xx[0] - roundToInt(0.5*(double)d[0]),
      y = xx[1] - roundToInt(0.5*(double)d[1]);

  XFillArc(XtDisplay(unixGUI.topLevel), unixGUI.pixmap, unixGUI.gc, 
	       x, y, d[0], d[1], a1, a2);
		  
//  XFillArc(XtDisplay(unixGUI.topLevel), XtWindow(unixGUI.drawingArea), unixGUI.gc,  
//	       x, y, d[0], d[1], a1, a2);
	
  return;
}




