
#include "Macro.h"
#include "DomainTree.h"


extern DomainTree domain;


int macro15(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "updt";
    macro.type = "anly";
    macro.what = "update domain for next time step";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;

    macro.db.selectDomain();

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int type, id;
  
  type = roundToInt(macro.p[0]);
  id   = roundToInt(macro.p[1]) - 1;

  domain(type,id).timeUpdate();

//--------------------------------------------------------------------------------------------------
  return 0;  
}

