
#include "Macro.h"
#include "DomainTree.h"
//#include "Plot.h"
//#include "Fluid.h"
//#include "Solid.h"


extern DomainTree domain;
//extern Plot       plot;


int macro49(Macro &macro)
{
/*
  if (!macro) 
  { 
    macro.name = "tran";
    macro.type = "anly";
    macro.what = "transfer data across interface";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;
    
    macro.db.selectDomain();

    macro.db.addRadioBox("*get displacments","transfer forces");

    macro.db.addRadioBox("*no relax.","relax. 1","relax. 2","relax. 3");

    macro.db.frameRadioBox();
    
    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int    type, id, what, how;
  
  type = roundToInt(macro.p[0]);
  id   = roundToInt(macro.p[1]) - 1;
  what = roundToInt(macro.p[2]);
  how  = roundToInt(macro.p[3]);
  
  if (what == 1) domain(type,id).getDisplacements();
 
  else           domain(type,id).transferForces();
*/
//--------------------------------------------------------------------------------------------------
  return 0;  
}
