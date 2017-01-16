
#include <Xm/Xm.h>
#include <Xm/Form.h>
#include <Xm/DialogS.h>
#include <Xm/PushB.h>
#include <Xm/LabelG.h>
#include <Xm/SeparatoG.h>

#include "RunControl.h"
#include "PlotVTK.h"

//extern UnixGUI unixGUI;
extern RunControl runCtrl;
extern PlotVTK  plotvtk;


//void grpReallyExit(Widget w, XtPointer client_data, XtPointer call_data);
//void grpDontExit  (Widget w, XtPointer client_data, XtPointer call_data);

	
//  exit mpap2

void grpExit(Widget w, XtPointer client_data, XtPointer call_data)
{
  runCtrl.quit = true;
  
  plotvtk.reset();
  
/*
  Widget dialogShell, baseForm, contrForm, separator, noButton, yesButton, message;

  XmString text = XmStringCreateSimple("Really exit mpap2 ?");

  dialogShell = XmCreateDialogShell(unixGUI.topLevel, "dialogShell", NULL, 0);  

  XtVaSetValues(dialogShell, XmNtitle, "MPAP2", myNULL);

  baseForm = XmCreateForm(dialogShell, "baseForm", NULL, 0);
  
  XtVaSetValues(baseForm, XmNdialogStyle, XmDIALOG_FULL_APPLICATION_MODAL, myNULL);
  
  message = XmCreateLabelGadget(baseForm, "message", NULL, 0);

  XtVaSetValues(message, XmNlabelString, text, 
		         XmNmarginBottom, 20,
		         XmNmarginTop,    20,
		         XmNmarginLeft,   20,
		         XmNmarginRight,  20,
			 XmNalignment,    XmALIGNMENT_CENTER,
			 //XmNfontList,    , 
		         //XmNwidth, 500, 
			 NULL);
  XmStringFree(text);
  
  separator = XmCreateSeparatorGadget(baseForm, "separator", NULL, 0); 

  contrForm = XmCreateForm(baseForm, "contrForm", NULL, 0);
  
  XtVaSetValues(message, XmNtopAttachment,    XmATTACH_FORM,
		         XmNbottomAttachment, XmATTACH_WIDGET,
			 XmNbottomWidget,     separator,
		         XmNleftAttachment,   XmATTACH_FORM,
		         XmNrightAttachment,  XmATTACH_FORM,
			 NULL);

  XtVaSetValues(separator, XmNbottomAttachment, XmATTACH_WIDGET,
			   XmNbottomWidget,     contrForm,
		           XmNleftAttachment,   XmATTACH_FORM,
		           XmNrightAttachment,  XmATTACH_FORM,
			   NULL);

  XtVaSetValues(contrForm, XmNbottomAttachment, XmATTACH_FORM,
		           XmNleftAttachment,   XmATTACH_FORM,
		           XmNrightAttachment,  XmATTACH_FORM,
		           XmNfractionBase,     35, 
			   NULL);
  
  yesButton = XmCreatePushButton(contrForm, "yes", NULL, 0);

  noButton = XmCreatePushButton(contrForm, "no", NULL, 0);

  XtVaSetValues(yesButton, XmNtopAttachment,                XmATTACH_POSITION,
	                   XmNtopPosition,                  5,
	                   XmNbottomAttachment,             XmATTACH_POSITION,
	                   XmNbottomPosition,               30,
	                   XmNleftAttachment,               XmATTACH_POSITION,   
	                   XmNleftPosition,                 7,
	                   XmNrightAttachment,              XmATTACH_POSITION,
	                   XmNrightPosition,                14,
	                   XmNshowAsDefault,                True,
	                   XmNdefaultButtonShadowThickness, 1,
		  	   NULL);
  
  XtVaSetValues(noButton, XmNtopAttachment,                XmATTACH_POSITION,
	                  XmNtopPosition,                  5,
	                  XmNbottomAttachment,             XmATTACH_POSITION,
	                  XmNbottomPosition,               30,
	                  XmNleftAttachment,               XmATTACH_POSITION,   
	                  XmNleftPosition,                 21,
	                  XmNrightAttachment,              XmATTACH_POSITION,
	                  XmNrightPosition,                28,
	                  XmNdefaultButtonShadowThickness, 1,
		  	  NULL);
  
  XtAddCallback(yesButton, XmNactivateCallback, grpReallyExit, baseForm);

  XtAddCallback(noButton, XmNactivateCallback, grpDontExit, baseForm);

  XtManageChild(yesButton);
  XtManageChild(noButton);
  XtManageChild(message);
  XtManageChild(separator);
  XtManageChild(contrForm);
  XtManageChild(baseForm);

}


void grpReallyExit(Widget w, XtPointer client_data, XtPointer call_data)
{
  grpDontExit(w, client_data, call_data);
      	
  exit_pressed = true;
}


void grpDontExit(Widget w, XtPointer client_data, XtPointer call_data)
{
  Widget baseForm = (Widget) client_data;
	
  XtUnmanageChild(baseForm);

  XtDestroyWidget(XtParent(baseForm));
*/
}


