
#include <cstdio>
#include <Xm/Xm.h>
#include <Xm/LabelG.h>

#include "FunctionsEssGrp.h"
#include "MpapTime.h"
#include "UnixGUI.h"


extern UnixGUI  unixGUI;
extern MpapTime mpapTime;
extern bool     noGUI;


void essGrpWriteTime(void)
{ 
  if (noGUI) return;

  Widget label = XtNameToWidget(unixGUI.topLevel,"*.timeLabel2");

  char tmp[40];

  XmString XmStr;

  sprintf(tmp, "%-12.5g", mpapTime.cur);
  
  int i = strlen(tmp) - 1;
  while (tmp[i] == ' ') tmp[i--] = '\0';

  XmStr = XmStringCreateSimple(tmp);

  XtVaSetValues(label, XmNlabelString, XmStr, myNULL);
  
  essGrpUpdateDisplay();
  
  XmStringFree(XmStr);

  return;
}


