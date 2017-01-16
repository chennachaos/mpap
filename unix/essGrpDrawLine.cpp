
#include <Xm/Xm.h>
#include <iostream>


#include "UnixGUI.h"


extern UnixGUI unixGUI;
extern bool    noGUI;


void essGrpDrawLine(int x11, int x12, int x21, int x22)
{
  if (noGUI) return;
      
  XDrawLine(XtDisplay(unixGUI.topLevel), unixGUI.pixmap, unixGUI.gc, x11, x12, x21, x22);
		  
//  XDrawLine(XtDisplay(unixGUI.topLevel), XtWindow(unixGUI.drawingArea), unixGUI.gc,  
//		  x11, x12, x21, x22);
	
  return;
}



