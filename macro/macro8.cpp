
#include "Macro.h"
#include "DomainTree.h"

extern DomainTree domain;


int macro8(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "disp";
    macro.type = "outp";
    macro.what = "print nodal displacements etc";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;
    macro.sensitivity[PRE]   = true;
    
    macro.db.selectDomain();
    
    macro.db.addTextField("node number = ",0,7);

    macro.db.addLabel("to print all nodal data, set to 0");
    

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int type, id, node;
  
  type = roundToInt(macro.p[0]);
  id   = roundToInt(macro.p[1]) - 1;
    
  node = roundToInt(macro.p[2]);
  
  domain(type,id).printNodalData(node);  

//--------------------------------------------------------------------------------------------------
  return 0;  
}

