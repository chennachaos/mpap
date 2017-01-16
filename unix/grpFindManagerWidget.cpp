
#include <Xm/Xm.h>

#include "FunctionsProgram.h"
#include "UnixGUI.h"


extern UnixGUI unixGUI;


Widget grpFindManagerWidget(Widget w)
{
  Widget      tmpWidget = w;

  while (strcmp(XtName(tmpWidget),"manager") != 0) 
    { 
       tmpWidget = XtParent(tmpWidget);
       if (tmpWidget == unixGUI.topLevel) prgError(1,"grpPopdownDialogShell","no 'manager' found!");
    }

  return tmpWidget;  
}


