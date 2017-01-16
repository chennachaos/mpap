
#include "Macro.h"
#include "Mesh.h"
#include "DomainTree.h"


extern DomainTree domain;


int macro74(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "sfix";
    macro.type = "prep";
    macro.what = "fix nodes";

    macro.sensitivity[PRE] = true;

    macro.db.selectDomain();

    macro.db.addToggleButton("fix x",0);
    macro.db.addToggleButton("fix y",0);
    macro.db.addToggleButton("fix z",0);

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int typ, id;

  bool fix[3];
 
  typ = roundToInt(macro.p[0]);
  id  = roundToInt(macro.p[1]) - 1;

  for (int i=0; i<3; i++) fix[i] = (roundToInt(macro.p[2+i]) == 1);

  mesh(domain(typ,id)).fixNodeMotion(fix);

//--------------------------------------------------------------------------------------------------
  return 0;  
}

