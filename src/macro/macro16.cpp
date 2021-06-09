
#include "Macro.h"
#include "DomainTree.h"


extern DomainTree domain;


int macro16(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "wrnd";
    macro.type = "outp";
    macro.what = "write nodal data to Tfile";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;

    macro.db.selectDomain();

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

    int i, type, id;

    type = roundToInt(macro.p[0]); 
    id   = roundToInt(macro.p[1]) - 1;

    domain(type,id).writeNodalData(); 

//--------------------------------------------------------------------------------------------------
  return 0;  
}

