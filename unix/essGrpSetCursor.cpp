
#include <iostream>
#include <Xm/Xm.h>
#include <X11/cursorfont.h>

#include "UnixGUI.h"
#include "RunControl.h"


extern UnixGUI    unixGUI;
extern bool       noGUI;
extern RunControl runCtrl;


void essGrpSetCursor(bool busy)
{
  if (noGUI) return;
  
  XSetWindowAttributes attrs;
  
  static Cursor busyCursor, pressMouseCursor, selectPosCursor, changePerspectiveCursor;

  if (busy)
  { 
    if (!busyCursor)  busyCursor = XCreateFontCursor(XtDisplay(unixGUI.topLevel), XC_X_cursor);
  
    attrs.cursor = busyCursor;
  }
  else   
  {
    if (runCtrl.fixedStatus == SELECTZOOMBOX
     || runCtrl.fixedStatus == SELECTNODE
     || runCtrl.fixedStatus == EXECTESTMACRO) 
    {
      if (!selectPosCursor) 
	selectPosCursor = XCreateFontCursor(XtDisplay(unixGUI.topLevel), XC_crosshair);
      
      attrs.cursor = selectPosCursor;
    }
    else if (runCtrl.fixedStatus == PRESSMOUSE)
    {
      if (!pressMouseCursor) 
	pressMouseCursor = XCreateFontCursor(XtDisplay(unixGUI.topLevel), XC_dot);
      
      attrs.cursor = pressMouseCursor;
    }
    else if (runCtrl.fixedStatus == CHANGEPERSPECTIVE)
    {
      if (!changePerspectiveCursor) 
	changePerspectiveCursor = XCreateFontCursor(XtDisplay(unixGUI.topLevel), XC_fleur);
      
      attrs.cursor = changePerspectiveCursor;
    }
    else 
    {
      attrs.cursor = None;
    }
  }
  
  XChangeWindowAttributes(XtDisplay(unixGUI.topLevel),XtWindow(unixGUI.topLevel), CWCursor, &attrs);

  return;
}



