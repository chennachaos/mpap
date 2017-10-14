
#include "Macro.h"
#include "DomainTree.h"


extern DomainTree domain;


int macro40(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "rset";
    macro.type = "anly";
    macro.what = "reset domain data to previous time instant";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;

    macro.db.selectDomain();

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int type = roundToInt(macro.p[0]),
      id   = roundToInt(macro.p[1]) - 1;
  
  domain(type,id).reset();
  
//--------------------------------------------------------------------------------------------------
  return 0;  
}

