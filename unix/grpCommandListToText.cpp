
#include <Xm/Xm.h>

#include "UnixGUI.h"



void grpCommandListToText(Widget w, XtPointer client_data, XtPointer call_data)
{
  char *tmp;
      
  Widget list = w, text = (Widget) client_data;

  XmListCallbackStruct *cbs = (XmListCallbackStruct*) call_data;
    
  XmStringGetLtoR(cbs->item, XmFONTLIST_DEFAULT_TAG, &tmp);   
  
  XtVaSetValues(text, XmNvalue, tmp, myNULL);
	
  delete [] tmp;
  
  return;
}



