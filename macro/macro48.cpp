
#include "Macro.h"
#include "DomainTree.h"
#include "Plot.h"
#include "Fluid.h"
#include "Solid.h"


extern DomainTree domain;
extern Plot       plot;


int macro48(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "intf";
    macro.type = "plot";
    macro.what = "plot nodes on interface boundary";

    macro.sensitivity[PRE] = true;
    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;
    
    macro.db.selectDomain();

    macro.db.addList(COLOURS_GREEN);

    macro.db.frameButtonBox();
    macro.db.frameRadioBox();

    macro.db.addRadioBox("free","*boundary","layers");

    macro.db.addRadioBox("initial","*current");

    macro.db.addToggleButton("node numbers", true);

    macro.db.addTextField("layers for node",0);
    
    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int    type, id, col, num, ndtyp, fnd;
  double xmn[3], xmx[3];
  bool   def;
  
  type  = roundToInt(macro.p[0]);
  id    = roundToInt(macro.p[1]) - 1;
  col   = roundToInt(macro.p[2]);
  ndtyp = roundToInt(macro.p[3]);
  def   =(roundToInt(macro.p[4]) == 2);
  num   = roundToInt(macro.p[5]);
  fnd   = roundToInt(macro.p[6]);
  
  plot.setColour(col-1);
  
  if (!plot) 
  {
     domain(type,id).findMinMaxX(xmn,xmx,def);

     plot.fit(xmn,xmx,domain(type,id).ndm);
  }

  domain(type,id).plotInterfaceNodes(num,ndtyp,fnd,def);

//--------------------------------------------------------------------------------------------------
  return 0;  
}

