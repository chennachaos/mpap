
#include <Xm/Xm.h>
#include <iostream>


#include "Debug.h"
#include "UnixGUI.h"


extern UnixGUI unixGUI;


void grpDrawingAreaExposed(Widget w, XtPointer client_data, XtPointer call_data)
{
  Dimension width, height;
  
  if (debug) std::cout << " DrawingArea has been exposed!\n\n";
  
  XtVaGetValues(unixGUI.drawingArea, XmNwidth, &width, XmNheight, &height, myNULL);
  
  XCopyArea(XtDisplay(unixGUI.topLevel), unixGUI.pixmap, XtWindow(unixGUI.drawingArea), 
		  unixGUI.gc,   0, 0, width, height, 0, 0);  

  return;
}

