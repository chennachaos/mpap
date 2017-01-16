
#include <Xm/Xm.h>
#include <Xm/Form.h>
#include <Xm/Frame.h>
#include <Xm/DialogS.h>
#include <Xm/PushB.h>
#include <Xm/LabelG.h>
#include <Xm/SeparatoG.h>

#include "Files.h"
#include "FunctionsProgram.h"
#include "FunctionsUnix.h"
#include "UnixGUI.h"


extern UnixGUI unixGUI;
extern Files  files;



void grpPopupLastProject(Widget w, XtPointer client_data, XtPointer call_data)
{
  Widget   dialogShell, 
	   manager, 
	   separator, 
	   frame, 
	   okButton, 
	   browseButton, 
	   cancelButton, 
	   message, 
	   project;
  
  XmString XmStr;
 
  char     tmp[500];
  
  dialogShell = XtVaCreatePopupShell("dialogShell", xmDialogShellWidgetClass, unixGUI.topLevel, 
		                     XmNtitle, "MPAP2", myNULL);
 
  manager    = XtVaCreateWidget("manager", xmFormWidgetClass, dialogShell,
		                XmNdialogStyle, XmDIALOG_FULL_APPLICATION_MODAL, myNULL);

  separator  = grpGenThreeButtonControlArea(&okButton,&browseButton,&cancelButton,manager,
	   	                            "OK","Browse","Cancel");

  XmStr = XmStringCreateLtoR("Run this project ?", XmSTRING_DEFAULT_CHARSET);
  
  message    = XtVaCreateManagedWidget("message", xmLabelGadgetClass, manager,
		                       XmNlabelString,      XmStr, 
   		                       XmNmarginBottom,     20,
		                       XmNmarginTop,        20,
		                       XmNmarginLeft,       20,
		                       XmNmarginRight,      20,
		                       XmNalignment,        XmALIGNMENT_CENTER,
                                       XmNbottomAttachment, XmATTACH_WIDGET,
			               XmNbottomWidget,     separator,
		                       XmNleftAttachment,   XmATTACH_FORM,
		                       XmNrightAttachment,  XmATTACH_FORM, myNULL);
  XmStringFree(XmStr);
  
  frame      = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, manager, 
		                       XmNleftOffset,       10,
		                       XmNrightOffset,      10,
		                       XmNtopOffset,        10, 
				       XmNtopAttachment,    XmATTACH_FORM,
		                       XmNbottomAttachment, XmATTACH_WIDGET,
			               XmNbottomWidget,     message,
		                       XmNleftAttachment,   XmATTACH_FORM,
		                       XmNrightAttachment,  XmATTACH_FORM, myNULL);

  if (files.projDir.length()
    + files.Ifile.length()
    + files.Ofile.length()
    + files.Tfile.length()
    + files.Pfile.length() > 400) 
	prgError(1,"grpLastProject","increase size of char tmp[500]!");

  sprintf(tmp,"directory: %s\n\n    Ifile: %s\n\n    Ofile: %s\n\n    Tfile: %s\n\n    Pfile: %s",
          files.projDir.asCharArray(), 
	  files.Ifile.asCharArray(), 
	  files.Ofile.asCharArray(), 
	  files.Tfile.asCharArray(), 
	  files.Pfile.asCharArray());
 
  XmStr = XmStringCreateLtoR(tmp, XmSTRING_DEFAULT_CHARSET);

  project     = XtVaCreateManagedWidget("project", xmLabelGadgetClass, frame, 
		                        XmNlabelString,  XmStr, 
  		                        XmNmarginBottom, 20,
		                        XmNmarginTop,    20,
		                        XmNmarginLeft,   20,
		                        XmNmarginRight,  20,
		                        XmNalignment,    XmALIGNMENT_BEGINNING, myNULL);
  XmStringFree(XmStr);

  XtAddCallback(okButton,     XmNactivateCallback, grpPopdownLastProject, 0);

  XtAddCallback(browseButton, XmNactivateCallback, grpPopdownLastProject, 0);
  
  XtAddCallback(cancelButton, XmNactivateCallback, grpPopdownLastProject, 0);

  XtManageChild(manager);
  
  return;
}





