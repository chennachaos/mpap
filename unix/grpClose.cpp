
#include <Xm/Xm.h>

#include "FunctionsEssGrp.h"
#include "FunctionsProgram.h"


void grpClose(Widget w, XtPointer client_data, XtPointer call_data)
{
  essGrpWipe();
	
  prgCloseProject();
	
  return;
}


