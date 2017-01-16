
#include "Macro.h"
#include "Mesh.h"
#include "DomainTree.h"


extern DomainTree domain;


int macro75(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "rfix";
    macro.type = "prep";
    macro.what = "free node motion";

    macro.sensitivity[PRE] = true;

    macro.db.selectDomain();

    macro.db.addToggleButton("free x",0);
    macro.db.addToggleButton("free y",0);
    macro.db.addToggleButton("free z",0);

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int typ, id;

  bool free[3];
 
  typ = roundToInt(macro.p[0]);
  id  = roundToInt(macro.p[1]) - 1;

  for (int i=0; i<3; i++) free[i] = (roundToInt(macro.p[2+i]) == 1);

  mesh(domain(typ,id)).fixNodeMotion(free,true);

//--------------------------------------------------------------------------------------------------
  return 0;  
}

