
#include <Xm/Xm.h>

#include "UnixGUI.h"
#include "RunControl.h"
#include "Plot.h"
#include "FunctionsEssGrp.h"
#include "FunctionsProgram.h"
#include "SelectNode.h"
#include "Domain.h"
#include "Perspective.h"


extern RunControl runCtrl;
extern UnixGUI    unixGUI;
extern Plot       plot;
extern void *macro2mousePtr;


using namespace std;


void grpDrawingAreaMouseInput(Widget w, XtPointer client_data, XtPointer call_data)
{
  Position x, y;
  XmDrawingAreaCallbackStruct *cbs = (XmDrawingAreaCallbackStruct*) call_data;
  XEvent *event = cbs->event;
  int dx, dy;
  double s1, s2;
  SelectNode *selectNode;

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// select zoom box
  
  if (runCtrl.fixedStatus == SELECTZOOMBOX)
  {
    if (event->xany.type == ButtonPress)
    { 
      x = event->xbutton.x;
      y = event->xbutton.y;
      plot.setColour(7);

      if (!plot.selectBox1) 
      {
        plot.selectBox1 = true;
        essGrpDrawLine(0,y,plot.wPix,y);
        essGrpDrawLine(x,0,x,plot.hPix);
	essGrpCopyPixmap();
        plot.sBx = x;
        plot.sBy = y;
      }      
      else if (x != plot.sBx && y != plot.sBy)
      {
        plot.selectBox1 = false;
        essGrpDrawLine(0,y,plot.wPix,y);
        essGrpDrawLine(x,0,x,plot.hPix);
        essGrpCopyPixmap();
        
	dx = abs(x-plot.sBx);
	dy = abs(y-plot.sBy);
	if (plot.sBx < x) x = plot.sBx;
	if (plot.sBy > y) y = plot.sBy;
      
        plot.x0Des[0] = ((double)x)             / plot.wFactAct + plot.x0Act[0];
        plot.x0Des[1] = ((double)(plot.hPix-y)) / plot.hFactAct + plot.x0Act[1];

        plot.dDes[0]  = ((double)dx) / plot.wFactAct;
        plot.dDes[1]  = ((double)dy) / plot.hFactAct;

        plot.adjustToNewSize();
	
	essGrpSetSensAllButDrawingArea(true);
	
	runCtrl.freeFixedStatus();
	runCtrl.notBusy();
      }
    }
    return;	  
  }

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// change 3D perspective (rotate, translate, zoom)

  Perspective &pers = plot.perspective;
  
  if (runCtrl.fixedStatus == CHANGEPERSPECTIVE)
  {
    if (event->xany.type == ButtonPress)
    { 
      //cout << " button " << event->xbutton.button << " pressed\n";

      plot.sBx = event->xbutton.x;
      plot.sBy = event->xbutton.y;
    }
    else if (event->xany.type == ButtonRelease)
    {
      //cout << " button " << event->xbutton.button << " released\n";

      if (event->xbutton.button == 2)
      {
        if (pers.setDesEqAct) 
        {
          plot.x0Des[0] = plot.x0Act[0];  plot.x0Des[1]  = plot.x0Act[1]; 
          plot.dDes[0]  = plot.dAct[0];   plot.dDes[1]   = plot.dAct[1];  
        }

        pers.print((float*)plot.dDes);

	essGrpSetSensAllButDrawingArea(true);
	
	runCtrl.freeFixedStatus();
	runCtrl.notBusy();
      }
      else
      {
        dx = event->xbutton.x - plot.sBx;
        dy = event->xbutton.y - plot.sBy;

        s1 = ((double) dx) / plot.wPix;
        s2 = ((double) dy) / plot.hPix;

        pers.calcNew((float)s1,(float)s2,plot.dAct,event->xbutton.button);

        essGrpWipe();
        plot.setColour(3);

        pers.plotChangePersObjSurf(pers.changePersObjSurfPtr->showBackFaces,
                                   pers.changePersObjSurfPtr->showDeformed);
	essGrpCopyPixmap();
      }
    }
    return;
  }
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// execute test macro
  
  if (runCtrl.fixedStatus == EXECTESTMACRO)
  {
    if (event->xany.type == ButtonPress)
    { 
      //cout << " button " << event->xbutton.button << " pressed\n";

      prgMacroTest((double) event->xbutton.x,(double) event->xbutton.y);

      essGrpSetSensAllButDrawingArea(true);
      runCtrl.freeFixedStatus();
      runCtrl.notBusy();
    }
    return;
  }
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// select node
  
  if (runCtrl.fixedStatus == SELECTNODE)
  {
    if (event->xany.type == ButtonPress)
    { 
      if (event->xbutton.button == 2)

        ((Domain*)macro2mousePtr)->interactiveNodeSelection(-1,-1);
 
      else

        ((Domain*)macro2mousePtr)->interactiveNodeSelection(
            (int)event->xbutton.x,(int)event->xbutton.y,event->xbutton.button == 3);

      if (event->xbutton.button == 2)
      {
        essGrpSetSensAllButDrawingArea(true);
        runCtrl.freeFixedStatus();
        runCtrl.notBusy();
      }
    }
    return;
  }
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  return;
}



