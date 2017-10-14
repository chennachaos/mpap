
#include "Macro.h"
#include "DomainTree.h"


extern DomainTree domain;


int macro29(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "proj";
    macro.type = "outp";
    macro.what = "project Gauss point data to nodes";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;
    macro.sensitivity[PRE]   = true;

    macro.db.selectDomain();

    macro.db.addTextField("index i",1);
    macro.db.addTextField("index j",1);

    macro.db.stringList("*stress",
                        "intvar",
                        "vorticity",
                        "error",
                        "gradient",
                        "normOfGradient",
                        "normOfGradientSquared"); 
    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int type, id, indx, indx2;
  
  type  = roundToInt(macro.p[0]);
  id    = roundToInt(macro.p[1]) - 1;
  indx  = roundToInt(macro.p[2]);
  indx2 = roundToInt(macro.p[2]);
 
  domain(type,id).projectToNodes(macro.strg,indx,indx2); 
 
//--------------------------------------------------------------------------------------------------
  return 0;  
}

