
#include <Xm/Xm.h>

#include "UnixGUI.h"
//#include "UnixGlobal.h"


extern UnixGUI unixGUI;
extern bool    noGUI;


void essGrpCopyPixmap(void)
{
  if (noGUI) return;

  Dimension width, height;
  
  XtVaGetValues(unixGUI.drawingArea, XmNwidth, &width, XmNheight, &height, myNULL);
  
  XCopyArea(XtDisplay(unixGUI.topLevel), unixGUI.pixmap, XtWindow(unixGUI.drawingArea), 
            unixGUI.gc,   0, 0, width, height, 0, 0);  

  return;
}



