
#include <Xm/Xm.h>
#include <Xm/DialogS.h>
#include <Xm/Form.h>
#include <Xm/PushB.h>
#include <Xm/LabelG.h>
#include <Xm/SeparatoG.h>
#include <Xm/RowColumn.h>
#include <Xm/ToggleB.h>
#include <Xm/List.h>
#include <Xm/TextF.h>
#include <Xm/BulletinB.h>
#include <Xm/Frame.h>
#include <cstdio>

#include "FunctionsUnix.h"
#include "FunctionsProgram.h"
#include "MacroList.h"
#include "MyString.h"
#include "DomainTree.h"
#include "MathBasic.h"
#include "UnixGUI.h"
#include "Definitions.h"


#define MAX_WIDGET 20



extern DomainTree domain;
extern MacroList  macro;
extern UnixGUI    unixGUI;



//  popup a macro dialog for input of macro parameters


void grpPopupMacro(Widget w, XtPointer client_data, XtPointer call_data)
{
  int      c, i, j, n, nargs, maci = (int) ((VOID_PTR) client_data), *ndom, ntype, dflt, 
	   cDflt = 0, cTxtF = 0, cTBtn = 0, cList = 0, 
	   cRBox = 0, cLabl = 0, cBttn = 0, cBBox = 0;

  Widget   dialogShell, 
	   manager,
           separator,
	   controlArea,
	   okButton, 
	   cancelButton, 
	   helpButton,
	   workArea,
	   selectDomain,
	   domainLabel,
	   macroLabel,
	   domainType,
	   domainNumber,
	   separator2,
	   boardWA,
	   boardSD,
	   strgList,
	   frame,
	   buttonBox,
	   list[MAX_WIDGET],
           rBox[MAX_WIDGET],
	   labl[MAX_WIDGET],
	   tBtn[MAX_WIDGET],
	   txtF[MAX_WIDGET],
	   txtFRC,
	   txtFLbl,
	   strgTxtF,
	   areaSTF,
	   boardSTF;

  MacroDialogData &db = macro[maci].db;

  XmString XmStr, *XmStrArray;
 
  char *tmp, tmp2[20];
      
  bool on;
  
  Arg args[4]; 

  Dimension w1, w2;
  
  if (db.list.n > MAX_WIDGET || db.rBox.n > MAX_WIDGET || db.labl.n > MAX_WIDGET
   || db.tBtn.n > MAX_WIDGET || db.txtF.n > MAX_WIDGET) 

    prgError(1,"grpPopupMacro","increase MAX_WIDGET!");

  
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// popup shell, manager form, all boards, all areas, all separators
	  
  dialogShell  = XtVaCreatePopupShell("dialogShell", xmDialogShellWidgetClass, unixGUI.topLevel, 
		                      XmNtitle,    "MPAP2", myNULL);
 
  manager      = XtVaCreateWidget("manager", xmFormWidgetClass, dialogShell,
		                  XmNdialogStyle, XmDIALOG_FULL_APPLICATION_MODAL, myNULL);
  
  controlArea  = grpGenThreeButtonControlArea(&okButton, &cancelButton, &helpButton, manager,
	 	                              "OK", "Cancel", "Help");

  XtVaSetValues(manager, XmNdefaultButton, XtNameToWidget(manager,"*.OK"),
                         XmNinitialFocus,  XtNameToWidget(manager,"*.controlForm"), myNULL);
  
  macroLabel   = XtVaCreateManagedWidget("macroLabel",xmLabelGadgetClass, manager,
		                         XmNmarginLeft,      15,
		                         XmNmarginRight,     15,
		                         XmNleftAttachment,  XmATTACH_FORM, 
		                         XmNrightAttachment, XmATTACH_FORM, myNULL);
  if (db.selDm) 
  {
    boardSD    = XtVaCreateManagedWidget("boardSD", xmBulletinBoardWidgetClass, manager,
		                         XmNorientation,     XmHORIZONTAL,
		                         XmNtopAttachment,   XmATTACH_FORM, 
		                         XmNleftAttachment,  XmATTACH_FORM, 
		                         XmNrightAttachment, XmATTACH_FORM, myNULL);

    selectDomain = XtVaCreateManagedWidget("selectDomain", xmRowColumnWidgetClass, boardSD,
              		                 XmNentryAlignment, XmALIGNMENT_CENTER,
		                         XmNorientation,    XmHORIZONTAL,
                                         XmNspacing,        15,
					 XmNmarginWidth,    20,
		                         XmNtopOffset,      10, myNULL);
 
    separator  = XtVaCreateManagedWidget("separator", xmSeparatorGadgetClass, manager,
		                         XmNseparatorType,    XmDOUBLE_LINE,
		                         XmNtopOffset,        5,
		                         XmNtopAttachment,    XmATTACH_WIDGET, 
		                         XmNtopWidget,        boardSD, 
		                         XmNleftAttachment,   XmATTACH_FORM, 
		                         XmNrightAttachment,  XmATTACH_FORM, myNULL);
    
    XtVaSetValues(macroLabel, XmNtopAttachment, XmATTACH_WIDGET, 
		              XmNtopWidget,     separator, myNULL);
  }
  else
        XtVaSetValues(macroLabel, XmNtopAttachment, XmATTACH_FORM, myNULL);

  if (db.inputType.dim() == 0 && db.strgList.n == 0 && !db.strgTxtFLabl)
  { 
        XtVaSetValues(controlArea, XmNseparatorType, XmDOUBLE_LINE, myNULL);

	XtVaSetValues(macroLabel, XmNbottomAttachment, XmATTACH_WIDGET, 
			          XmNbottomWidget,     controlArea, myNULL);
  }
  else 
  {	  
    separator2 = XtVaCreateManagedWidget("separator2", xmSeparatorGadgetClass, manager,
		                         XmNseparatorType,   XmDOUBLE_LINE,
		                         XmNtopAttachment,   XmATTACH_WIDGET, 
		                         XmNtopWidget,       macroLabel, 
		                         XmNleftAttachment,  XmATTACH_FORM, 
		                         XmNrightAttachment, XmATTACH_FORM, myNULL);
  
    if (db.inputType.dim() > 0 || db.strgList.n > 0)  
    {
      boardWA  = XtVaCreateManagedWidget("boardWA", xmBulletinBoardWidgetClass, manager,
		                         XmNtopAttachment,     XmATTACH_WIDGET, 
		                         XmNtopWidget,         separator2, 
		                         XmNleftAttachment,    XmATTACH_FORM, 
		                         XmNrightAttachment,   XmATTACH_FORM, myNULL);
    
      workArea = XtVaCreateManagedWidget("workArea", xmRowColumnWidgetClass, boardWA,
		                         XmNentryAlignment, XmALIGNMENT_CENTER,
		                         XmNorientation,    XmHORIZONTAL,
		                         XmNspacing,        35, myNULL);
    }    
    if (!!db.strgTxtFLabl)
    {
      boardSTF = XtVaCreateManagedWidget("boardSTF", xmBulletinBoardWidgetClass, manager,
		                         XmNbottomAttachment,  XmATTACH_WIDGET, 
		                         XmNbottomWidget,      controlArea, 
		                         XmNleftAttachment,    XmATTACH_FORM, 
		                         XmNrightAttachment,   XmATTACH_FORM, myNULL);
      
      areaSTF  = XtVaCreateManagedWidget("areaSTF", xmRowColumnWidgetClass, boardSTF,
		                         XmNentryAlignment, XmALIGNMENT_CENTER,
		                         XmNorientation,    XmHORIZONTAL, myNULL);

      if (db.inputType.dim() > 0 || db.strgList.n > 0)
	      
	  XtVaSetValues(boardSTF, XmNtopAttachment, XmATTACH_WIDGET, 
			          XmNtopWidget,     boardWA, myNULL);
      else
	  XtVaSetValues(boardSTF, XmNtopAttachment, XmATTACH_WIDGET,
			          XmNtopWidget,     separator2, myNULL);
    }
    else 
          XtVaSetValues(boardWA, XmNbottomAttachment, XmATTACH_WIDGET,
		                 XmNbottomWidget,     controlArea, myNULL);
  }
  
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// manage domain selection and prepare macro label
 
  if (db.selDm) 
  {
    XmStr = XmStringCreateSimple("select domain:");
   
    domainLabel = XtVaCreateManagedWidget("domainLabel", xmLabelGadgetClass, selectDomain,
     	                                  XmNlabelString, XmStr, myNULL); 
    XmStringFree(XmStr);

    ntype      = domain.nType();
    XmStrArray = new XmString[ntype];
    ndom       = new int     [ntype];
    c = 0;
    for (i=0; i<=domain.mxtid; i++) 
      {  j = domain.nDomainOfType(i);
         if (j > 0) 
           {  ndom[c]       = j; 
     	      XmStrArray[c] = XmStringCreateSimple(domain[i].type.asCharArray());
     	      c++;
           }
      }
    nargs = 0;
    XtSetArg(args[nargs], XmNbuttons,     XmStrArray); nargs++;
    XtSetArg(args[nargs], XmNbuttonCount, ntype);      nargs++;
    XtSetArg(args[nargs], XmNbuttonSet,   0);          nargs++;
    
    domainType = XmCreateSimpleOptionMenu(selectDomain, "domainType", args, nargs);

    for (int i=0; i<ntype; i++) XmStringFree(XmStrArray[i]);
    delete [] XmStrArray;
    delete [] ndom;

    XtManageChild(domainType);
   
    domainNumber = XtVaCreateManagedWidget("domainNumber", xmTextFieldWidgetClass, selectDomain, 
	                                   XmNcolumns, 5, 
		                           XmNvalue, "1", myNULL);
  }

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// set macro label string

  tmp = new char[4+4+10+macro[maci].what.length()];
  sprintf(tmp,"%s,%s - %s",macro[maci].type.asCharArray(),
		           macro[maci].name.asCharArray(),
			   macro[maci].what.asCharArray());
  XmStr = XmStringCreateSimple(tmp);
  delete [] tmp;
  
  XtVaSetValues(macroLabel, XmNlabelString, XmStr, 
		            XmNheight,      50, myNULL);

  XmStringFree(XmStr);

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// string list or string text field

  if (db.strgList.n > 0)
    {
       n = db.strgList.n;
       
       XmStrArray = new XmString[n];

       for (j=0;j<n;j++) 
	 XmStrArray[j] = XmStringCreateSimple(db.strgList[j].asCharArray());

       strgList = XmCreateScrolledList(workArea, "strgList", NULL, 0);
       
       XtVaSetValues(strgList,XmNitems,            XmStrArray,
			      XmNitemCount,        n,
			      XmNvisibleItemCount, 4, myNULL);

       XmStr = XmStringCreateSimple(db.dfltStrg.asCharArray());
       XmListSelectItem(strgList,XmStr,false);
       XmStringFree(XmStr);

       XtManageChild(strgList);
	    
       for (j=0;j<n;j++) XmStringFree(XmStrArray[j]);
       
       delete [] XmStrArray; 
    } 

  if (!!db.strgTxtFLabl)
    {
      XmStr = XmStringCreateSimple(db.strgTxtFLabl.asCharArray());

      txtFLbl = XtVaCreateManagedWidget("txtFLbl", xmLabelGadgetClass, areaSTF,
               				XmNmarginLeft,  5,
               				XmNlabelString, XmStr, myNULL);
      XmStringFree(XmStr);  
      
      strgTxtF = XtVaCreateManagedWidget("strgTxtF", xmTextFieldWidgetClass, areaSTF, 
               				             XmNvalue,   db.dfltStrg.asCharArray(),
               				             XmNcolumns, db.strgTxtFCol, myNULL); 
    }
  
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//

  //if (debug) for (i=0; i<db.nBttn.dim(); i++) std::cout << " nBttn " << db.nBttn[i] << "\n\n";
  
  cBBox = -1;
  
  tmp = new char[20];
 
  for (i=0; i<db.inputType.dim(); i++)
  {
    // generate new button box if necessary

    if (db.inputType[i] != LIST && db.inputType[i] != RBOX)

      if (cBttn == 0 || cBttn == db.nBttn[cBBox])
      {
	 cBBox++; cBttn = 0;

         c = 1; while (c*db.bttnCol[cBBox]<db.nBttn[cBBox]) c++;
 
         if (db.frmBBoxFlg)
	  
         {  frame = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, workArea, myNULL);
	                     
            buttonBox = XtVaCreateManagedWidget("ButtonBox", xmRowColumnWidgetClass, frame,
		                                XmNrowColumnType, XmWORK_AREA,
				                XmNpacking,       XmPACK_COLUMN,
				                XmNnumColumns,    c,
		                                XmNorientation,   XmVERTICAL, myNULL);   }
         else
            buttonBox = XtVaCreateManagedWidget("ButtonBox", xmRowColumnWidgetClass, workArea,
		                                XmNrowColumnType, XmWORK_AREA,
				                XmNpacking,       XmPACK_COLUMN,
				                XmNnumColumns,    c,
		                                XmNorientation,   XmVERTICAL, myNULL);	 
      } 

    // generate input widget
    
    switch (db.inputType[i])
    {
      case LIST: // generate list
	      
	         //std::cout << " List " << cList << ",   InputType (" << i << ")\n\n";
		 
                 n = db.list[cList].n;
                 
                 XmStrArray = new XmString[n];

                 for (j=0;j<n;j++) 
               	   XmStrArray[j] = XmStringCreateSimple(db.list[cList][j].asCharArray());

                 sprintf(tmp,"list%d",cList);
                 
                 list[cList] = XmCreateScrolledList(workArea, tmp, NULL, 0);
                 
                 XtVaSetValues(list[cList],XmNitems,            XmStrArray,
               	   	                   XmNitemCount,        n,
               	   	                   XmNvisibleItemCount, 4, myNULL);
                
		 dflt = roundToInt(db.dflt[cDflt++]);

                 XmListSelectPos(list[cList], dflt, false);

                 XtManageChild(list[cList]);
               	 
                 for (j=0;j<n;j++) XmStringFree(XmStrArray[j]);
                 
                 delete [] XmStrArray; 

		 cList++;

		 break;

      case RBOX: // generate radio box
		 
		 //std::cout << " Radio Box " << cRBox << ",   InputType (" << i << ")\n\n";
		 
                 n = db.rBox[cRBox].n;
                      
                 XmStrArray = new XmString[n];

                 for (j=0;j<n;j++) 
               	   XmStrArray[j] = XmStringCreateSimple(db.rBox[cRBox][j].asCharArray());

		 dflt = roundToInt(db.dflt[cDflt++]) - 1;

                 nargs = 0;
                 XtSetArg(args[nargs], XmNbuttons,     XmStrArray); nargs++;
                 XtSetArg(args[nargs], XmNbuttonCount, n);          nargs++;
                 XtSetArg(args[nargs], XmNbuttonSet,   dflt);       nargs++;
                 
                 sprintf(tmp,"rBox%d",cRBox);

                 if (db.frmRBoxFlg)
            	     
                 {  frame = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, workArea, myNULL);
            	       
                    rBox[cRBox] = XmCreateSimpleRadioBox(frame, tmp, args, nargs);   }
            	 
                 else
                    rBox[cRBox] = XmCreateSimpleRadioBox(workArea, tmp, args, nargs);

                 XtManageChild(rBox[cRBox]);
            	  
                 for (j=0;j<n;j++) XmStringFree(XmStrArray[j]);
                 
                 delete [] XmStrArray; 

		 cRBox++;

		 break;

      case LABL: // generate label in current button box
		
                 XmStr = XmStringCreateSimple(db.labl[cLabl].asCharArray());

                 sprintf(tmp,"labl%d",cLabl);

               	 labl[cLabl] = XtVaCreateManagedWidget(tmp, xmLabelGadgetClass, buttonBox, 
               			                            XmNlabelString, XmStr, myNULL);
                 XmStringFree(XmStr);  

                 cLabl++; cBttn++;
		 
		 break;

      case TBTN: // generate toggle button in current button box

                 XmStr = XmStringCreateSimple(db.tBtn[cTBtn].asCharArray());

                 sprintf(tmp,"tBtn%d",cTBtn);

                 on = false; if (roundToInt(db.dflt[cDflt++]) == 1) on = true;
		 
                 tBtn[cTBtn] = XtVaCreateManagedWidget(tmp, xmToggleButtonWidgetClass, buttonBox, 
               				                    XmNset,         on,
               			                            XmNlabelString, XmStr, myNULL);
                 XmStringFree(XmStr);  

                 cTBtn++; cBttn++;
		 
		 break;

      case TXTF: // generate text field in current button box

                 XmStr = XmStringCreateSimple(db.txtF[cTxtF].asCharArray());

                 sprintf(tmp,"txtF%d",cTxtF);

                 txtFRC  = XtVaCreateManagedWidget("txtFRC", xmRowColumnWidgetClass, buttonBox,
               		                                     XmNorientation, XmHORIZONTAL, myNULL);
               		     
    		 txtFLbl = XtVaCreateManagedWidget("txtFLbl", xmLabelGadgetClass, txtFRC,
               				                      XmNmarginLeft,  5,
               				                      XmNlabelString, XmStr, myNULL);
               		    
                 sprintf(tmp2,"%g",db.dflt[cDflt++]); 

                 txtF[cTxtF] = XtVaCreateManagedWidget(tmp, xmTextFieldWidgetClass, txtFRC, 
               				                    XmNvalue,   tmp2,
               				                    XmNcolumns, db.txtFCol[cTxtF], myNULL); 
                 XmStringFree(XmStr);  

                 cTxtF++; cBttn++;
		 
		 break;
    }
  }
  delete [] tmp;
  
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// add callbacks and manage everything

  XtAddCallback(okButton,     XmNactivateCallback, grpPopdownMacro, (XtPointer) maci);
  
  XtAddCallback(cancelButton, XmNactivateCallback, grpPopdownMacro, (XtPointer) maci);
  
  XtAddCallback(helpButton,   XmNactivateCallback, grpMacroHelp, (XtPointer) macro[maci].id);

  XtManageChild(manager);

  if (db.selDm) 
    {
       XtVaGetValues(selectDomain, XmNwidth, &w1, myNULL);
       XtVaGetValues(boardSD,      XmNwidth, &w2, myNULL);

       XtVaSetValues(selectDomain, XmNx, (w2-w1)/2, myNULL);
    }
  
  if (db.inputType.dim() > 0 || db.strgList.n > 0) 
    {  
       XtVaGetValues(workArea, XmNwidth, &w1, myNULL);
       XtVaGetValues(boardWA,  XmNwidth, &w2, myNULL);

       XtVaSetValues(workArea, XmNx, (w2-w1)/2, myNULL);
    }

  if (!!db.strgTxtFLabl)
    {
       XtVaGetValues(areaSTF,  XmNwidth, &w1, myNULL);
       XtVaGetValues(boardSTF, XmNwidth, &w2, myNULL);

       XtVaSetValues(areaSTF,  XmNx, (w2-w1)/2, myNULL);
    }
 
  return;
}


