
#include <iostream>

#include <Xm/Xm.h>

#include "FunctionsProgram.h"
#include "Definitions.h"


extern bool noGUI;


void grpMacroHelp(Widget w, XtPointer client_data, XtPointer call_data)
{
  if (noGUI) return;
	
  MyString helpCmd;

  char tmp[10];

  int macid = (int) ((VOID_PTR) client_data); 
 
  sprintf(tmp,"%d",macid);

  helpCmd.free().append("mpap2help macros")
                .append(SLASH)
                .append("macro")
                .append(tmp)
                .append(".html");
  
  system(helpCmd.asCharArray());
	  
  return;
}
