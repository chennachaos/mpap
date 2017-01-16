
#include "Macro.h"
#include "DomainTree.h"


extern DomainTree domain;



int macro99(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "lift";
    macro.type = "wulf";
    macro.what = "lifting line analysis of finite wing";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;

    macro.db.selectDomain();

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int type, id;
  
  type = roundToInt(macro.p[0]);
  id   = roundToInt(macro.p[1]) - 1;

  domain(type,id).doForLiftingLine();

//--------------------------------------------------------------------------------------------------
  return 0;  
}

