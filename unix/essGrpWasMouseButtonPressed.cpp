
#include <Xm/Xm.h>

#include "UnixGUI.h"
#include "FunctionsEssGrp.h"
#include "FunctionsProgram.h"


extern UnixGUI unixGUI;
extern bool    noGUI;


Bool mouseButtonPressed(Display*, XEvent*, char*);


int essGrpWasMouseButtonPressed(int *x, int *y, bool *inDrawingArea)
{
  if (noGUI) return 0;

  XEvent event;
 
  int key = 0;
 
  if (XCheckIfEvent(XtDisplay(unixGUI.topLevel), &event, mouseButtonPressed, NULL))
  {
    if (event.xany.type == ButtonPress)   
    {
      key = (int) event.xbutton.button;

      if (x != NULL) *x = (int) event.xbutton.x;
      if (y != NULL) *y = (int) event.xbutton.y;

      if (inDrawingArea != NULL) 
        *inDrawingArea = (event.xbutton.window == XtWindow(unixGUI.drawingArea));
    }
    
    else prgError(1,"essGrpWasMouseButtonPressed","event type mismatch!");
  }
  
  if (event.xany.type == ConfigureNotify) 
    if (event.xany.window == XtWindow(unixGUI.topLevel)) XtDispatchEvent(&event);

  return key;
}





Bool mouseButtonPressed(Display*, XEvent *event, char*)
{
  if (event->xany.type == ButtonPress) return true;

  return false;
}


