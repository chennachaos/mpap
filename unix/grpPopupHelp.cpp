
#include <Xm/Xm.h>

#include "FunctionsProgram.h"


extern bool noGUI;


void grpPopupHelp(Widget w, XtPointer client_data, XtPointer call_data)
{
  if (noGUI) return;
	
  system("mpap2help content.html");
	  
  return;
}
