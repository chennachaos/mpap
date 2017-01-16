
#include <Xm/Xm.h>

#include "MacroList.h"
#include "RunControl.h"
#include "UnixGUI.h"


extern MacroList  macro;
extern bool       noGUI;
extern UnixGUI    unixGUI;
extern RunControl runCtrl;


void essGrpSetMacroSens(void)
{
  if (noGUI) return;
  
  char tmp[20];
  
/*  if (runCtrl.mode == NOPROJ)
  {	  
    for (int i=0; i<macro.n; i++)
    {    
       sprintf(tmp,"*.%s",macro[i].name.asCharArray());
  	  
       XtVaSetValues(XtNameToWidget(unixGUI.topLevel,tmp), XmNsensitive, false, myNULL);
    }
    return;
  }
*/
  for (int i=0; i<macro.n; i++)
  {    
     sprintf(tmp,"*.%s",macro[i].name.asCharArray());
	  
     XtVaSetValues(XtNameToWidget(unixGUI.topLevel,tmp), XmNsensitive, 
		     macro[i].sensitivity[runCtrl.mode], myNULL);
  }

  return;
}



