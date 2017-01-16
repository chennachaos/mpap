
#include <Xm/Xm.h>
#include <Xm/Form.h>
#include <Xm/PushB.h>
#include <Xm/LabelG.h>
#include <Xm/SeparatoG.h>

#include "UnixGUI.h"


Widget grpGenTwoButtonControlArea(Widget *PBW1, Widget *PBW2, Widget managerForm,
		                  char   *PBL1, char   *PBL2)
{
  Widget controlForm, separator;
  
  controlForm  = XtVaCreateManagedWidget("controlForm", xmFormWidgetClass, managerForm,
                                         XmNbottomAttachment, XmATTACH_FORM,
	                                 XmNleftAttachment,   XmATTACH_FORM,
	                                 XmNrightAttachment,  XmATTACH_FORM,  
	                                 XmNfractionBase,     27, myNULL);
 
  separator    = XtVaCreateManagedWidget("separator2", xmSeparatorGadgetClass, managerForm,
       		                         XmNbottomAttachment, XmATTACH_WIDGET,
	       	                         XmNbottomWidget,     controlForm,
	                                 XmNleftAttachment,   XmATTACH_FORM,
	                                 XmNrightAttachment,  XmATTACH_FORM, myNULL);

  *PBW1        = XtVaCreateManagedWidget(PBL1, xmPushButtonWidgetClass, controlForm,
	         			 XmNshadowThickness,  3,
	                                 XmNmarginHeight,     6,
                                         XmNtopAttachment,    XmATTACH_POSITION, 
	        	                 XmNtopPosition,      5,
	                                 XmNbottomAttachment, XmATTACH_POSITION, 
	        	                 XmNbottomPosition,   21,
	                                 XmNleftAttachment,   XmATTACH_POSITION, 
	        	                 XmNleftPosition,     5,
	                                 XmNrightAttachment,  XmATTACH_POSITION, 
	        	                 XmNrightPosition,    11, myNULL);
          
  *PBW2        = XtVaCreateManagedWidget(PBL2, xmPushButtonWidgetClass, controlForm,
	         			 XmNshadowThickness,  3,
                                         XmNtopAttachment,    XmATTACH_POSITION, 
	        	                 XmNtopPosition,      5,
	                                 XmNbottomAttachment, XmATTACH_POSITION, 
	        	                 XmNbottomPosition,   21,
	                                 XmNleftAttachment,   XmATTACH_POSITION, 
	        	                 XmNleftPosition,     16,
	                                 XmNrightAttachment,  XmATTACH_POSITION, 
	        	                 XmNrightPosition,    22, myNULL);
  return separator;
}


