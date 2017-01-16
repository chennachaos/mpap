
#include <Xm/Xm.h>

#include "Definitions.h"
#include "RunControl.h"
#include "UnixGUI.h"


extern RunControl runCtrl;
extern UnixGUI    unixGUI;
extern bool       noGUI;


void essGrpWriteStatus(void)
{ 
  if (noGUI) return;

  Widget label = XtNameToWidget(unixGUI.topLevel,"*.statusLabel");
  
  char *statusName[] = STATUS_NAMES;
      
  XmString XmStr;

  if (runCtrl.fixedStatus != UNDEF) XmStr = XmStringCreateSimple(statusName[runCtrl.fixedStatus]);
  else                              XmStr = XmStringCreateSimple(statusName[runCtrl.status]);
  
  XtUnmanageChild(label);

  XtVaSetValues(label, XmNlabelString, XmStr, myNULL);
  
  XtManageChild(label);

  essGrpUpdateDisplay();
  
  XmStringFree(XmStr);

  return;
}


