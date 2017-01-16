
#include <Xm/Xm.h>

#include "UnixGUI.h"
#include "FunctionsEssGrp.h"

extern UnixGUI unixGUI;
extern bool    noGUI;


int essGrpWaitForMouseButtonPressed(void)
{
  if (noGUI) return 0;

  XEvent event;
 
  int key = 0;
 
  essGrpSetCursor(false);
  
  while (key < 1)
  {
    XNextEvent(XtDisplay(unixGUI.topLevel), &event);
    
    if (event.xany.type == ButtonPress) key = (int) event.xbutton.button;

    if (event.xany.type == ConfigureNotify) 
      if (event.xany.window == XtWindow(unixGUI.topLevel)) XtDispatchEvent(&event);
  }

  essGrpSetCursor(true);

  return key;
}



