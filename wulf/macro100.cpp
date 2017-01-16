
#include "Macro.h"
#include "DomainTree.h"
#include "Mesh.h"


extern DomainTree domain;


int macro100(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "vsm";
    macro.type = "wulf";
    macro.what = "vortex sheet method in 2D";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;

    macro.db.selectDomain();

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int type, id;
  
  type = roundToInt(macro.p[0]);
  id   = roundToInt(macro.p[1]) - 1;

  domain(type,id).doForVortexSheet();

//--------------------------------------------------------------------------------------------------
  return 0;  
}

