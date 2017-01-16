
#include <iostream>
#include <Xm/Xm.h>

#include "UnixGUI.h"
//#include "UnixGlobal.h"


extern UnixGUI unixGUI;
extern bool    noGUI;


void essGrpCheckForResizeEvents(void)
{
  if (noGUI) return;

  XEvent event;
  
  while (XCheckMaskEvent(XtDisplay(unixGUI.topLevel), StructureNotifyMask, &event)) 
  {
    if (event.xany.window == XtWindow(unixGUI.topLevel))  
    {
      //std::cout << " StructureNotifyMask:  event dispatched!\n\n";	    
      
      XtDispatchEvent(&event);
    }
  }
  
  return;
}






