
#include <Xm/Xm.h>
#include <Xm/DialogS.h>
#include <Xm/FileSB.h>
#include <Xm/LabelG.h>
#include <Xm/Form.h>
#include <Xm/SeparatoG.h>
#include <Xm/Text.h>

#include "FunctionsProgram.h"
#include "FunctionsUnix.h"
#include "Files.h"
#include "Definitions.h"
#include "MyString.h"
#include "Debug.h"
#include "UnixGUI.h"


void grpPopdownOutFileNames(Widget, XtPointer, XtPointer);

	
extern Files files;
extern UnixGUI unixGUI;


void grpPopdownFileSelect(Widget w, XtPointer client_data, XtPointer call_data)
{
  Widget fileNameText, argWidget = (Widget) client_data;

  if (strcmp(XtName(argWidget),"Cancel") == 0) 
  
         { grpPopdownDialogShell(grpFindManagerWidget(w)); return; }


  // retrieve projDir and Ifile from fileSelection

  char *tmp;

  int  l;
  
  fileNameText = XmFileSelectionBoxGetChild(argWidget, XmDIALOG_TEXT);
 
  if (!(tmp = XmTextGetString(fileNameText))) return;

  if (!prgFileExist(tmp)) return;

  l = strlen(tmp)-1; while (tmp[l] != SLASH) l--;
       
  files.projDir.free().append(tmp).trunc(l);
       
  files.Ifile.free().append(&(tmp[l+1]));

  XtFree(tmp);
 
  grpPopdownDialogShell(grpFindManagerWidget(w));

  
  // popup dialog to input output file names

  tmp = new char[files.Ifile.length() + 3];
  
  Widget dialogShell, 
	 manager, 
	 label, 
	 separator, 
	 okButton, 
	 cancelButton, 
	 controlArea, 
	 textOfile,
	 textTfile, 
	 textPfile,
	 labelOfile,
	 labelTfile, 
	 labelPfile;
 
  XmString XmStr;
  
  dialogShell   = XtVaCreatePopupShell("dialogShell", xmDialogShellWidgetClass, unixGUI.topLevel, 
		                       XmNtitle, "MPAP2", myNULL);

  manager       = XtVaCreateWidget("manager", xmFormWidgetClass, dialogShell,
		                   XmNverticalSpacing, 10,
		                   XmNdialogStyle,  XmDIALOG_FULL_APPLICATION_MODAL, myNULL);
  
  XmStr = XmStringCreateLtoR("Input output file names!", XmSTRING_DEFAULT_CHARSET);

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

  controlArea   = grpGenTwoButtonControlArea(&okButton, &cancelButton, manager, "OK", "Cancel");

  XmStr = XmStringCreateLtoR("output file:", XmSTRING_DEFAULT_CHARSET);
  
  labelOfile    = XtVaCreateManagedWidget("labelOfile", xmLabelGadgetClass, manager,
		                          XmNlabelString,      XmStr,
					  XmNwidth,            180,
					  XmNmarginLeft,       20,
					  XmNheight,           35,
					  XmNalignment,        XmALIGNMENT_BEGINNING,
		                          XmNtopAttachment,    XmATTACH_WIDGET,
					  XmNtopWidget,        separator,
					  XmNleftAttachment,   XmATTACH_FORM, myNULL);
  XmStringFree(XmStr);	  	  
  
  sprintf(tmp,"O%s",&(files.Ifile.asCharArray()[1]));
  
  textOfile     = XtVaCreateManagedWidget("textOfile", xmTextWidgetClass, manager,
					  XmNvalue,            tmp,
					  XmNheight,           35,
					  XmNrightOffset,      10,
		                          XmNtopAttachment,    XmATTACH_WIDGET,
					  XmNtopWidget,        separator,
					  XmNleftAttachment,   XmATTACH_WIDGET,
					  XmNleftWidget,       labelOfile,
					  XmNrightAttachment,  XmATTACH_FORM, myNULL);

  XmStr = XmStringCreateLtoR("time plot file:", XmSTRING_DEFAULT_CHARSET);
  
  labelTfile    = XtVaCreateManagedWidget("labelTfile", xmLabelGadgetClass, manager,
					  XmNheight,           35,
					  XmNwidth,            180,
					  XmNmarginLeft,       20,
					  XmNalignment,        XmALIGNMENT_BEGINNING,
		                          XmNlabelString,      XmStr,
		                          XmNtopAttachment,    XmATTACH_WIDGET,
					  XmNtopWidget,        labelOfile,
					  XmNleftAttachment,   XmATTACH_FORM, myNULL);
  XmStringFree(XmStr);	  	  
  
  sprintf(tmp,"T%s",&(files.Ifile.asCharArray()[1]));
  
  textTfile     = XtVaCreateManagedWidget("textTfile", xmTextWidgetClass, manager,
					  XmNvalue,            tmp,
					  XmNrightOffset,      10,
					  XmNheight,           35,
		                          XmNtopAttachment,    XmATTACH_WIDGET,
					  XmNtopWidget,        textOfile,
					  XmNleftAttachment,   XmATTACH_WIDGET,
					  XmNleftWidget,       labelTfile,
					  XmNrightAttachment,  XmATTACH_FORM, myNULL);
  
  XmStr = XmStringCreateLtoR("eps plot file:", XmSTRING_DEFAULT_CHARSET);
  
  labelPfile    = XtVaCreateManagedWidget("labelPfile", xmLabelGadgetClass, manager,
					  XmNheight,           35,
					  XmNwidth,            180,
					  XmNmarginLeft,       20,
		                          XmNlabelString,      XmStr,
					  XmNalignment,        XmALIGNMENT_BEGINNING,
		                          XmNtopAttachment,    XmATTACH_WIDGET,
					  XmNtopWidget,        labelTfile,
					  XmNleftAttachment,   XmATTACH_FORM, myNULL);
  XmStringFree(XmStr);	  	  
  
  sprintf(tmp,"P%s",&(files.Ifile.asCharArray()[1]));

  textPfile     = XtVaCreateManagedWidget("textPfile", xmTextWidgetClass, manager,
					  XmNrightOffset,      10,
					  XmNheight,           35,
					  XmNvalue,            tmp,
		                          XmNtopAttachment,    XmATTACH_WIDGET,
					  XmNtopWidget,        textTfile,
		                          XmNbottomAttachment, XmATTACH_WIDGET,
					  XmNbottomWidget,     controlArea,
					  XmNleftAttachment,   XmATTACH_WIDGET,
					  XmNleftWidget,       labelPfile,
					  XmNrightAttachment,  XmATTACH_FORM, myNULL);
  delete [] tmp;

  XtAddCallback(okButton,     XmNactivateCallback, grpPopdownOutFileNames, manager);
  
  XtAddCallback(cancelButton, XmNactivateCallback, grpPopdownOutFileNames, 0);
  
  XtManageChild(manager);
  
  return;  
}





void grpPopdownOutFileNames(Widget w, XtPointer client_data, XtPointer call_data)
{
  if (strcmp(XtName(w),"Cancel") == 0) 

    { grpPopdownDialogShell(grpFindManagerWidget(w)); 
	    
      files.Ifile.free(); 
      files.projDir.free(); 
      
      return; 
    }

  char *tmp;
  
  Widget text, manager = (Widget) client_data;
  
  text = XtNameToWidget(manager, "*.textOfile");
  if (!(tmp = XmTextGetString(text))) return;
  files.Ofile = tmp;  
  XtFree(tmp);
  if (!files.Ofile.FileOK()) return;

  text = XtNameToWidget(manager, "*.textTfile");
  if (!(tmp = XmTextGetString(text))) return;
  files.Tfile = tmp;  
  XtFree(tmp);
  if (!files.Tfile.FileOK()) return;

  text = XtNameToWidget(manager, "*.textPfile");
  if (!(tmp = XmTextGetString(text))) return;
  files.Pfile = tmp;  
  XtFree(tmp);
  if (!files.Pfile.FileOK()) return;

  grpPopdownDialogShell(grpFindManagerWidget(w));

  if (debug) 
    {  std::cout << " " << files.projDir << "\n";
       std::cout << " " << files.Ifile << "\n";
       std::cout << " " << files.Ofile << "\n";
       std::cout << " " << files.Tfile << "\n";
       std::cout << " " << files.Pfile << "\n\n"; }
  
  prgWriteLastProject();

  grpOpenProject();
  
  return;
}






