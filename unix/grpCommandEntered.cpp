
#include <Xm/Xm.h>
#include <Xm/List.h>
#include <iostream>


#include "MacroCommand.h"
#include "FunctionsProgram.h"
#include "MyString.h"
#include "MacroQueue.h"
#include "UnixGUI.h"


extern MacroQueue macroQueue;


void grpCommandEntered(Widget w, XtPointer client_data, XtPointer call_data)
{
  Widget list, text;
  
  XmListCallbackStruct *cbs;
  
  char *tmp;  
  
  int  nCmd;
  
  MyString myStr;
  
  XmString XmStr;

  MacroCommand macCmd;
  
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// retrieve the command string from list or text widget
  
  if (strcmp(XtName(w),"commandText") == 0)
  {
    list = (Widget) client_data;
    text = w;
    
    XtVaGetValues(text, XmNvalue, &tmp, myNULL);
  }
  else
  {
    list = w;
    text = (Widget) client_data;
    
    cbs = (XmListCallbackStruct*) call_data;
    
    XmStringGetLtoR(cbs->item, XmFONTLIST_DEFAULT_TAG, &tmp);   
  }
  
  myStr = tmp;

  delete [] tmp;

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// generate macro command from input string, return if not admissible

  if (!prgStringToMacroCmd(macCmd,myStr)) return;

  XtVaSetValues(text, XmNvalue, "", myNULL);
  
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// improve the format of the macro command and add it to the list
 
  prgMacroCmdToSimpleString(myStr,macCmd);

  XmStr = XmStringCreateSimple(myStr.asCharArray());
  
  XmListAddItemUnselected(list,XmStr,0);
  XmListSetBottomPos(list,0);
	
  XmStringFree(XmStr); 

  XtUnmanageChild(XtParent(list));
  XtManageChild(XtParent(list));
  
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// add the macro to the macro command queue 

  macroQueue.append(macCmd);       

  return;
}

