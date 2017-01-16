
#include <Xm/Xm.h>
#include <Xm/List.h>
#include <Xm/ToggleB.h>
#include <Xm/Command.h>

#include "FunctionsProgram.h"
#include "FunctionsUnix.h"
#include "MacroCommand.h"
#include "MacroList.h"
#include "MacroDialogData.h"
#include "DomainTree.h"
#include "DomainTypeEnum.h"
#include "MacroQueue.h"
#include "UnixGUI.h"
#include "Definitions.h"


extern DomainTree domain;
extern MacroList  macro;
extern UnixGUI    unixGUI;
extern MacroQueue macroQueue;


void grpPopdownMacro(Widget w, XtPointer client_data, XtPointer call_data)
{
	
  Widget dmy, manager = grpFindManagerWidget(w);

  int maci = (int) ((VOID_PTR) client_data); 
  
  if (strcmp(XtName(w),"Cancel") == 0) { grpPopdownDialogShell(manager); return; }

  char *tmp, *tmp2;

  int i, j, ip, c, pc = 0, *selPos, pos, cList = 0, cRBox = 0, cTBtn = 0, cTxtF = 0;

  double dp;
  
  bool on;

  MacroCommand macCmd;

  MyString macCmdStr;
 
  XmString XmStr;
  
  macCmd.ii = maci;

  MacroDialogData &db = macro[maci].db;

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// read domain selection
  
  if (db.selDm) 
  {
    XtVaGetValues(XtNameToWidget(manager, "*.domainType"), XmNmenuHistory, &dmy, myNULL);  

    tmp = new char[20]; 
    
    sprintf(tmp,XtName(dmy));
    
    scanInt(&tmp[7],&ip);
    
    delete [] tmp;

    c = 0; 
    i = 0;
    while (1) { if (domain.nDomainOfType(i) > 0) { if (c == ip) break; c++; } i++; }
    
    macCmd.p[pc++].x = (double) i;
    
    XtVaGetValues(XtNameToWidget(manager, "*.domainNumber"), XmNvalue, &tmp, myNULL);

    if (!macCmd.p[pc++].interprete(tmp,NON_NEG_INT)) { delete [] tmp; return; }
    
    if (roundToInt(macCmd.p[pc-1].x) > domain[i].dom.n)  macCmd.p[pc-1].x = domain[i].dom.n;
    
    delete [] tmp;
  }

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// read string list or string text field 

  if (db.strgList.n > 0) 
  {
    XmListGetSelectedPos(XtNameToWidget(manager, "*.strgList"), &selPos, &c);
    pos = selPos[0] - 1;
    delete [] selPos;
      
    macCmd.strg.free().append(db.strgList[pos]);
  }

  if (!!db.strgTxtFLabl)
  {
    XtVaGetValues(XtNameToWidget(manager, "*.strgTxtF"), XmNvalue, &tmp, myNULL);
	 
    j = 0; while (!charInCharArray(tmp[j], ALL_WORD_SEPARATORS)) j++; 
    if (j<strlen(tmp)) { delete [] tmp; return; }
    
    macCmd.strg.free().append(tmp); 
    
    delete [] tmp;
  }

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// read input devices

  tmp = new char[20];
 
  for (i=0; i<db.inputType.dim(); i++)
  {
    switch (db.inputType[i])
    {
      case LIST: sprintf(tmp,"*.list%d",cList++);
      
                 XmListGetSelectedPos(XtNameToWidget(manager, tmp), &selPos, &c);

                 macCmd.p[pc++].x = (double) selPos[0];

                 delete [] selPos;

		 break;

      case RBOX: sprintf(tmp,"*.rBox%d",cRBox++);
      
                 dmy = XtNameToWidget(manager,tmp);

                 pos = 0;
                 while (1)
                 {
                    sprintf(tmp,"*.button_%d",pos);
	   
                    XtVaGetValues(XtNameToWidget(dmy,tmp), XmNset, &on, myNULL);
      
                    pos++; if (on) break; 
                 }
      
                 macCmd.p[pc++].x = (double) pos;
		 
		 break;

      case TBTN: sprintf(tmp,"*.tBtn%d",cTBtn++);
		 
		 XtVaGetValues(XtNameToWidget(manager,tmp), XmNset, &on, myNULL);
	       
                 if (on) macCmd.p[pc++].x = 1;
		 else    macCmd.p[pc++].x = 0;

		 break;

      case TXTF: sprintf(tmp,"*.txtF%d",cTxtF++);
		 
		 XtVaGetValues(XtNameToWidget(manager,tmp), XmNvalue, &tmp2, myNULL);

                 if (!macCmd.p[pc++].interprete(tmp2,ALL_REAL)) { delete [] tmp2; return; }
		 
                 delete [] tmp2;

		 break;
    }
  }
  delete [] tmp;

  
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// prepare macro command string and pass it to the command widget

  prgMacroCmdToSimpleString(macCmdStr,macCmd);
  
  // add the macro with the parameters to the command widget
     
  dmy = XtNameToWidget(unixGUI.topLevel, "*.commandList");
  
  XmStr = XmStringCreateSimple(macCmdStr.asCharArray());

  XmListAddItemUnselected(dmy,XmStr,0); 

  XmStringFree(XmStr);

  XmListSetBottomPos(dmy,0);

  XtVaSetValues(XtNameToWidget(unixGUI.topLevel, "*.commandText"), XmNvalue, "", myNULL);

  XtUnmanageChild(XtParent(dmy));
  XtManageChild(XtParent(dmy));
  
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// popdown the macro dialog
  
  grpPopdownDialogShell(manager);

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// add the macro to the macro command queue 

  macroQueue.append(macCmd);       



  return;  
}


