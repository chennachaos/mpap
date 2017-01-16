
#include <Xm/Xm.h>
#include <iostream>


#include "UnixGUI.h"
#include "FunctionsEssGrp.h"


extern bool    noGUI;
extern UnixGUI unixGUI;


void essGrpWipe(void)
{
  if (noGUI) return;

  Dimension width, height;
      
  // wipe pixmap

  XtVaGetValues(unixGUI.drawingArea, XmNwidth, &width, XmNheight, &height, myNULL);
  
  XSetForeground(XtDisplay(unixGUI.drawingArea), unixGUI.gc, 
		 BlackPixelOfScreen (XtScreen(unixGUI.drawingArea)));

  XFillRectangle(XtDisplay(unixGUI.drawingArea), unixGUI.pixmap, unixGUI.gc, 0, 0, width, height);  

  return;
}
  
