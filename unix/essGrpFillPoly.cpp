
#include <Xm/Xm.h>
#include <iostream>


#include "UnixGUI.h"


extern UnixGUI unixGUI;
extern bool    noGUI;



void essGrpFillPoly(int *pnt, int npt)
{
  if (noGUI) return;
      
  if (npt > unixGUI.nXPoints)
  {
    delete [] unixGUI.XPoints; 
    unixGUI.XPoints  = new XPoint [npt]; 
    unixGUI.nXPoints = npt; 
  }
	
  for (int i=0; i<npt; i++)
  {
    unixGUI.XPoints[i].x = pnt[i+i];
    unixGUI.XPoints[i].y = pnt[i+i+1];
  }
	
  XFillPolygon(XtDisplay(unixGUI.topLevel), unixGUI.pixmap, unixGUI.gc, 
	       unixGUI.XPoints, npt, Complex, CoordModeOrigin);
		  
//  XFillPolygon(XtDisplay(unixGUI.topLevel), XtWindow(unixGUI.drawingArea), unixGUI.gc,  
//	       unixGUI.XPoints, npt, Complex, CoordModeOrigin);
	
  return;
}




