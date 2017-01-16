
#include <iostream>

#include <Xm/Xm.h>
#include <Xm/DialogS.h>
#include <Xm/Form.h>
#include <Xm/PushB.h>
#include <Xm/LabelG.h>
#include <Xm/SeparatoG.h>
#include <Xm/TextF.h>

#include "FunctionsUnix.h"
#include "UnixGUI.h"
#include "MyString.h"


void grpInquireButton1(Widget, XtPointer, XtPointer);
void grpInquireButton2(Widget, XtPointer, XtPointer);


extern UnixGUI unixGUI;


bool essGrpInquire(char* title, char* bttn1, char* bttn2, MyString* strg)
{
  int answer = 0;
	
  Widget inquireShell, 
         manager,
	 titleLabel,
         button1,
         button2,
	 textField,
	 topOfControlForm;
 
  XmString XmStr;
 
  char *tmp;
	 
  inquireShell = XtVaCreatePopupShell("inquireShell", xmDialogShellWidgetClass,
                                    unixGUI.topLevel, XmNtitle, "MPAP2", myNULL);
 
  manager      =   XtVaCreateWidget("manager", xmFormWidgetClass, inquireShell,
                                    XmNdialogStyle,  XmDIALOG_FULL_APPLICATION_MODAL, myNULL);
  
  topOfControlForm = grpGenTwoButtonControlArea(&button1, &button2, manager, bttn1, bttn2);
  
  XmStr = XmStringCreateSimple(title);
  
  titleLabel = XtVaCreateManagedWidget("titleInquire", xmLabelGadgetClass, manager,
      	                            XmNlabelString,     XmStr, 
				    XmNheight,          50,
      		  	            XmNtopAttachment,   XmATTACH_FORM,
      			            XmNleftAttachment,  XmATTACH_FORM,
      			            XmNrightAttachment, XmATTACH_FORM,
				    XmNleftOffset,      10,
				    XmNrightOffset,     10, myNULL);
  XmStringFree(XmStr);
  
  if (strg != NULL)
    textField = XtVaCreateManagedWidget("textField", xmTextFieldWidgetClass, manager, 
                                    XmNvalue,            strg->asCharArray(),
				    XmNleftAttachment,   XmATTACH_FORM,
				    XmNrightAttachment,  XmATTACH_FORM,
				    XmNtopAttachment,    XmATTACH_WIDGET,
				    XmNbottomAttachment, XmATTACH_WIDGET,
				    XmNtopWidget,        titleLabel,
				    XmNbottomWidget,     topOfControlForm, 
				    XmNleftOffset,       50,  
				    XmNrightOffset,      50,  
				    XmNbottomOffset,     10, myNULL); 
  else
    XtVaSetValues(topOfControlForm,XmNtopAttachment, XmATTACH_WIDGET,
		                   XmNtopWidget,     titleLabel, myNULL);
  
 
  XtAddCallback(button1,  XmNactivateCallback, grpInquireButton1, &answer);
  XtAddCallback(button2,  XmNactivateCallback, grpInquireButton2, &answer);
  
  XtManageChild(manager);

  while (answer == 0)  XtAppProcessEvent(unixGUI.app,XtIMAll);
  
  grpPopdownDialogShell(manager);

  if (answer == 1)
  {
    XtVaGetValues(textField, XmNvalue, &tmp, myNULL);

    strg->free();
    strg->append(tmp);

    delete [] tmp;
	  
    return true;
  }

  return false;
}




void grpInquireButton1(Widget w, XtPointer client_data, XtPointer call_data)
{
  int *answer = (int*) client_data;  *answer = 1;  return;
}



void grpInquireButton2(Widget w, XtPointer client_data, XtPointer call_data)
{
  int *answer = (int*) client_data;  *answer = 2;  return;
}


