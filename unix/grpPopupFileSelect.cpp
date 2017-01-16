
#include <Xm/Xm.h>
#include <Xm/DialogS.h>
#include <Xm/FileSB.h>
#include <Xm/Form.h>
#include <Xm/SeparatoG.h>
#include <Xm/LabelG.h>

#include "FunctionsProgram.h"
#include "FunctionsUnix.h"
#include "UnixGUI.h"


extern UnixGUI unixGUI;


void grpPopupFileSelect(Widget w, XtPointer client_data, XtPointer call_data)
{
  Widget dialogShell, fileSelection, manager, label, separator, button;
 
  XmString XmStr, XmStr2, XmStr3;
  
  dialogShell   = XtVaCreatePopupShell("dialogShell", xmDialogShellWidgetClass, unixGUI.topLevel, 
		                       XmNtitle, "MPAP2", myNULL);


  manager       = XtVaCreateWidget("manager", xmFormWidgetClass, dialogShell,
		                   XmNdialogStyle, XmDIALOG_FULL_APPLICATION_MODAL, myNULL);
  
  XmStr = XmStringCreateLtoR("Select the MPAP2 input file!", XmSTRING_DEFAULT_CHARSET);

  label         = XtVaCreateManagedWidget("label", xmLabelGadgetClass, manager,
		                          XmNmarginHeight,     20,
		                          XmNlabelString,      XmStr,
					  XmNtopAttachment,    XmATTACH_FORM,
					  XmNleftAttachment,   XmATTACH_FORM,
					  XmNrightAttachment,  XmATTACH_FORM, myNULL);
  XmStringFree(XmStr);
  
  separator     = XtVaCreateManagedWidget("separator", xmSeparatorGadgetClass, manager,
		                          XmNtopAttachment,    XmATTACH_WIDGET,
					  XmNtopWidget,        label, 
					  XmNleftAttachment,   XmATTACH_FORM,
					  XmNrightAttachment,  XmATTACH_FORM, myNULL);

  XmStr  = XmStringCreateSimple("~/Documents/IGACode/mpap3/project");
  XmStr2 = XmStringCreateSimple(" [  ...  ]          ");
  XmStr3 = XmStringCreateSimple("I*");
  
  fileSelection = XtVaCreateManagedWidget("fileSelect", xmFileSelectionBoxWidgetClass, manager,
		                      XmNdirectory,        XmStr,
                          XmNnoMatchString,    XmStr2,
                          XmNpattern,          XmStr3,
                          XmNtopAttachment,    XmATTACH_WIDGET,
                          XmNtopWidget,        separator,
                          XmNbottomAttachment, XmATTACH_FORM,
                          XmNleftAttachment,   XmATTACH_FORM,
                          XmNrightAttachment,  XmATTACH_FORM, 
                          XmNwidth,            480, myNULL);
  
  XmStringFree(XmStr);
  XmStringFree(XmStr2);
  XmStringFree(XmStr3);
  
  //XtUnmanageChild(XmFileSelectionBoxGetChild(fileSelection,XmDIALOG_DIR_LIST));
  XtUnmanageChild(XmFileSelectionBoxGetChild(fileSelection,XmDIALOG_HELP_BUTTON));

  button        = XmFileSelectionBoxGetChild(fileSelection,XmDIALOG_OK_BUTTON);
 
  XtVaSetValues(button, XmNdefaultButtonShadowThickness, 0, 
                        XmNshadowThickness,  3,
                        XmNmarginHeight,     6, myNULL);
  
  button        = XmFileSelectionBoxGetChild(fileSelection,XmDIALOG_APPLY_BUTTON);
 
  XtVaSetValues(button, XmNdefaultButtonShadowThickness, 0, 
                        XmNshadowThickness,  3,
                        XmNmarginHeight,     6, myNULL);
  
  button        = XmFileSelectionBoxGetChild(fileSelection,XmDIALOG_CANCEL_BUTTON);
 
  XtVaSetValues(button, XmNdefaultButtonShadowThickness, 0, 
                        XmNshadowThickness,  3,
                        XmNmarginHeight,     6, myNULL);
 
  XtAddCallback(fileSelection, XmNokCallback,     grpPopdownFileSelect, fileSelection);

  XtAddCallback(fileSelection, XmNcancelCallback, grpPopdownFileSelect, button);

  XtManageChild(manager);

  return;
}



