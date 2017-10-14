
#include "Macro.h"
//#include "Mesh.h"
#include "DomainTree.h"


extern DomainTree domain;


int macro67(Macro &macro)
{
/*
  if (!macro) 
  { 
    macro.name = "sbc";
    macro.type = "prep";
    macro.what = "set boundary conditions";

    macro.sensitivity[PRE] = true;

    macro.db.selectDomain();

    macro.db.addToggleButton("fix 1st d.o.f.",0);
    macro.db.addToggleButton("fix 2nd d.o.f.",0);
    macro.db.addToggleButton("fix 3rd d.o.f.",0);
    macro.db.addToggleButton("fix 4th d.o.f.",0);
    macro.db.addToggleButton("fix 5th d.o.f.",0);
    macro.db.addToggleButton("fix 6th d.o.f.",0);

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int typ, id;

  bool fix[20];
 
  typ = roundToInt(macro.p[0]);
  id  = roundToInt(macro.p[1]) - 1;

  for (int i=0; i<6; i++) fix[i] = (roundToInt(macro.p[2+i]) == 1);

  mesh(domain(typ,id)).setBoundaryConditions(fix);
*/
//--------------------------------------------------------------------------------------------------
  return 0;  
}

