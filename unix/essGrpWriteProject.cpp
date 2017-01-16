
#include <Xm/Xm.h>

#include "Files.h"
#include "UnixGUI.h"


extern UnixGUI unixGUI;
extern Files   files;
extern bool    noGUI;


void essGrpWriteProject(char *strg)
{ 
  if (noGUI) return;
	
  XmString XmStr;

  if (strg != NULL)
    {
       XmStr = XmStringCreateSimple(strg);
       XtVaSetValues(XtNameToWidget(unixGUI.topLevel,"*.open"),  XmNsensitive, false, myNULL);
       XtVaSetValues(XtNameToWidget(unixGUI.topLevel,"*.close"), XmNsensitive, true, myNULL);
    }
  else if (!files.Ifile)  
    {
       XmStr = XmStringCreateSimple(" ");
       XtVaSetValues(XtNameToWidget(unixGUI.topLevel,"*.open"),  XmNsensitive, true, myNULL);
       XtVaSetValues(XtNameToWidget(unixGUI.topLevel,"*.close"), XmNsensitive, false, myNULL);
    }
  else
    {
       XmStr = XmStringCreateSimple(files.Ifile.asCharArray());
       XtVaSetValues(XtNameToWidget(unixGUI.topLevel,"*.open"),  XmNsensitive, false, myNULL);
       XtVaSetValues(XtNameToWidget(unixGUI.topLevel,"*.close"), XmNsensitive, true, myNULL);
    }

  XtUnmanageChild(XtNameToWidget(unixGUI.topLevel,"*.projNameLabel"));

  XtVaSetValues(XtNameToWidget(unixGUI.topLevel,"*.projNameLabel"), XmNlabelString, XmStr, myNULL);
  
  XtManageChild(XtNameToWidget(unixGUI.topLevel,"*.projNameLabel"));

  XmStringFree(XmStr);

  return;
}


