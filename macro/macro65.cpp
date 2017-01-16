
#include "Macro.h"
#include "RunControl.h"
//#include "Mesh.h"
#include "DomainTree.h"


extern RunControl runCtrl;
extern DomainTree domain;
extern void *macro2mousePtr;


int macro65(Macro &macro)
{
/*
  if (!macro) 
  { 
    macro.name = "ins";
    macro.type = "prep";
    macro.what = "interactive node selection";

    macro.sensitivity[PRE]   = true;
    macro.sensitivity[INTER] = true;

    macro.db.selectDomain();

    macro.db.addRadioBox("*geometry mode","single nodes");

    macro.db.addTextField   ("max angle = ",20);
    macro.db.addToggleButton("print nodelist",1);
    macro.db.addToggleButton("show numbers",0);

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int typ, id;
 
  typ = roundToInt(macro.p[0]);
  id  = roundToInt(macro.p[1]) - 1;

  SelectNode &selectNode = mesh(domain(typ,id)).selectNode;

  selectNode.flag     = roundToInt(macro.p[2]);
  selectNode.alpha    = min(180.,macro.p[3]);
  selectNode.cosAlpha = cos(selectNode.alpha*0.017453293);
  selectNode.prnt     = (roundToInt(macro.p[4]) == 1);
  selectNode.numb     = (roundToInt(macro.p[5]) == 1);

  macro2mousePtr = (void*)&domain(typ,id);

  essGrpSetSensAllButDrawingArea(false);
   
  runCtrl.fixStatus(SELECTNODE);
*/
//--------------------------------------------------------------------------------------------------
  return 0;  
}

