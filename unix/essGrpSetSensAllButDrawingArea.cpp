
#include <Xm/Xm.h>

#include "UnixGUI.h"
#include "MacroList.h"


extern bool      noGUI;
extern UnixGUI   unixGUI;
extern MacroList macro;


void essGrpSetSensAllButDrawingArea(bool Sens)
{
  if (noGUI) return;

  char tmp[100];
 
  // macro menu
  
  for (int i=0; i<macro.ntype; i++)
  {    
    sprintf(tmp,"*.menuBar.%s",macro[macro.jp[i]].type.asCharArray());

    XtSetSensitive(XtNameToWidget(unixGUI.topLevel,tmp), Sens);
  }

  // project, help
  
  XtSetSensitive(XtNameToWidget(unixGUI.topLevel,"*.menuBar.project"), Sens);
  XtSetSensitive(XtNameToWidget(unixGUI.topLevel,"*.menuBar.help"),    Sens);

  // command list and text field

  XtSetSensitive(XtNameToWidget(unixGUI.topLevel,"*.commandList"), Sens);
  XtSetSensitive(XtNameToWidget(unixGUI.topLevel,"*.commandText"), Sens);
  
  return;
}





