
#include <Xm/Xm.h>
#include <iostream>


#include "Debug.h"
#include "UnixGUI.h"
#include "FunctionsEssGrp.h"
#include "Plot.h"


extern UnixGUI unixGUI;
extern bool    noGUI;
extern Plot    plot;


void grpDrawingAreaResized(Widget w, XtPointer client_data, XtPointer call_data)
{
  if (noGUI) return;
  
  Dimension width, height;
  
  if (debug) std::cout << " DrawingArea has been resized!\n\n";

  // free pixmap
  
  XFreePixmap(XtDisplay(unixGUI.topLevel), unixGUI.pixmap);
	
  // generate new pixmap
  
  XtVaGetValues(unixGUI.drawingArea, XmNwidth, &width, XmNheight, &height, myNULL);
  
  unixGUI.pixmap = XCreatePixmap(XtDisplay(unixGUI.drawingArea), 
		         RootWindowOfScreen(XtScreen(unixGUI.drawingArea)), width, height,
 		         DefaultDepthOfScreen(XtScreen(unixGUI.drawingArea)));

  // wipe the screen and the pixmap

  essGrpWipe();
 
  essGrpCopyPixmap();
  
  plot.adjustToNewSize();

  // make sure the scrolled list of commands looks ok
  
  Widget cmdListSW = XtNameToWidget(unixGUI.topLevel, "*.cmdListSW");
  XtUnmanageChild(cmdListSW);
  XtManageChild(cmdListSW);
  

  // after each 'resized' event there automatically follows an 'exposed' event.

  
  return;
}

