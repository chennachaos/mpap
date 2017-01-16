
#include <Xm/Xm.h>
#include <Xm/DialogS.h>
#include <Xm/Form.h>
#include <Xm/PushB.h>
#include <Xm/LabelG.h>
#include <Xm/SeparatoG.h>

#include "FunctionsUnix.h"
#include "UnixGUI.h"


extern UnixGUI unixGUI;

void grpPopdownAbout(Widget, XtPointer, XtPointer);


//#include "FunctionsProgram.h"


void grpPopupAbout(Widget w, XtPointer client_data, XtPointer call_data)
{
  Widget aboutShell,
	 manager,
	 closeButton,
	 label1, label2, separ;	
 
  XmString XmStr;
  
  XmFontList      fontList;
  XFontStruct     *fontStruct;
  XmFontListEntry fontListEntry;  
  
  aboutShell  = XtVaCreatePopupShell("aboutShell", xmDialogShellWidgetClass, unixGUI.topLevel, 
	 	                     XmNtitle, "MPAP2", myNULL);
 
  manager     = XtVaCreateWidget("manager", xmFormWidgetClass, aboutShell,
				 XmNmarginHeight,   20,
				 XmNmarginWidth,    10,
	 	                 XmNdialogStyle,    XmDIALOG_FULL_APPLICATION_MODAL, myNULL);
  
  fontStruct    = XLoadQueryFont(XtDisplay(unixGUI.topLevel),"*helvetica-bold-r-normal--34*");
  fontListEntry = XmFontListEntryCreate("tag1", XmFONT_IS_FONT, fontStruct);
  fontList      = XmFontListAppendEntry(NULL, fontListEntry);

  fontStruct    = XLoadQueryFont(XtDisplay(unixGUI.topLevel),"*helvetica-medium-r-normal--20*");
  fontListEntry = XmFontListEntryCreate("tag2", XmFONT_IS_FONT, fontStruct);
  fontList      = XmFontListAppendEntry(fontList, fontListEntry);

  XmStr  = XmStringCreate("M P A P 2","tag1");
		  
  label1 = XtVaCreateManagedWidget("label1", xmLabelGadgetClass, manager,
				   XmNleftAttachment, XmATTACH_FORM,
				   XmNrightAttachment, XmATTACH_FORM,
				   XmNtopOffset,     10,
		                   XmNfontList,      fontList,
		                   XmNlabelString,   XmStr, 
				   XmNtopAttachment, XmATTACH_FORM, myNULL);
  XmStringFree(XmStr);  
  
  XmStr  = XmStringCreate(" MULTI - PHYSICS  ANALYSIS  PROGRAM ","tag2");
		  
  label2 = XtVaCreateManagedWidget("label2", xmLabelGadgetClass, manager,
				   XmNleftAttachment, XmATTACH_FORM,
				   XmNrightAttachment, XmATTACH_FORM,
				   XmNtopOffset,     10,
		                   XmNfontList,      fontList,
		                   XmNlabelString,   XmStr,  
				   XmNtopAttachment, XmATTACH_WIDGET, 
				   XmNtopWidget,     label1, myNULL);
  XmStringFree(XmStr);  

  separ  = XtVaCreateManagedWidget("separator", xmSeparatorGadgetClass, manager,
				   XmNtopOffset,     5,
				   XmNseparatorType, XmSINGLE_LINE,
				   XmNleftAttachment, XmATTACH_FORM,
				   XmNrightAttachment, XmATTACH_FORM,
		                   XmNtopAttachment, XmATTACH_WIDGET, 
				   XmNtopWidget,     label2, myNULL);

  XmStr  = XmStringCreate(" University  of  Wales  Swansea ","tag2");
		  
  label1 = XtVaCreateManagedWidget("label2", xmLabelGadgetClass, manager,
				   XmNtopOffset,     5,
				   XmNleftAttachment, XmATTACH_FORM,
				   XmNrightAttachment, XmATTACH_FORM,
		                   XmNfontList,      fontList,
		                   XmNlabelString,   XmStr,  
				   XmNtopAttachment, XmATTACH_WIDGET, 
				   XmNtopWidget,     separ, myNULL);
  XmStringFree(XmStr);  

  XmStr  = XmStringCreate("WGD 2006","tag2");
		  
  label2 = XtVaCreateManagedWidget("label2", xmLabelGadgetClass, manager,
				   XmNtopOffset,     30,
				   XmNleftAttachment, XmATTACH_FORM,
				   XmNrightAttachment, XmATTACH_FORM,
				   XmNbottomAttachment, XmATTACH_FORM,
		                   XmNfontList,      fontList,
		                   XmNlabelString,   XmStr,  
				   XmNtopAttachment, XmATTACH_WIDGET, 
				   XmNtopWidget,     label1, myNULL);
  XmStringFree(XmStr);  
  
  XmFontListFree(fontList);
  
  XmStr = XmStringCreateSimple("close");
  
  closeButton = XtVaCreateManagedWidget("closeAbout", xmPushButtonWidgetClass, manager,
		                        XmNlabelString,   XmStr,  
				        XmNrightAttachment,  XmATTACH_FORM,
				        XmNbottomAttachment, XmATTACH_FORM,
	         		        XmNshadowThickness,  3,
                                        XmNmarginHeight,     6, myNULL);
          
  XmStringFree(XmStr);  
  
  XtAddCallback(closeButton, XmNactivateCallback, grpPopdownAbout, 0);

  XtManageChild(manager);
  
  return;
}





void grpPopdownAbout(Widget w, XtPointer client_data, XtPointer call_data)
{
  grpPopdownDialogShell(grpFindManagerWidget(w));

  return;
}


