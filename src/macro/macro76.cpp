
#include "Macro.h"
#include "DomainTree.h"
//#include "Mesh.h"


extern DomainTree domain;


int macro76(Macro &macro)
{
/*
  if (!macro) 
  { 
    macro.name = "spds";
    macro.type = "prep";
    macro.what = "set prescribed displacements";

    macro.sensitivity[PRE] = true;

    macro.db.selectDomain();
 
    macro.db.addList("*uniform","parabolic");

    macro.db.addToggleButton("1st d.o.f.",0);
    macro.db.addToggleButton("2nd d.o.f.",0);
    macro.db.addToggleButton("3rd d.o.f.",0);
    macro.db.addToggleButton("4th d.o.f.",0);
    macro.db.addToggleButton("5th d.o.f.",0);
    macro.db.addToggleButton("6th d.o.f.",0);

    macro.db.addTextField("",0,6);
    macro.db.addTextField("",0,6);
    macro.db.addTextField("",0,6);
    macro.db.addTextField("",0,6);
    macro.db.addTextField("",0,6);
    macro.db.addTextField("",0,6);

    macro.db.nextButtonBox();

    macro.db.addTextField("xp",0,6);
    macro.db.addTextField("yp",0,6);
    macro.db.addTextField("zp",0,6);
    macro.db.addTextField("r ",0,6);

    macro.db.setButtonColumnDim(6,4);

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int typ, id, flag;

  bool fix[6];

  double val[6], xp[3], rad;
 
  typ  = roundToInt(macro.p[0]);
  id   = roundToInt(macro.p[1]) - 1;
  flag = roundToInt(macro.p[2]);

  for (int i=0; i<6; i++)
  {
    fix[i] = (roundToInt(macro.p[3+i]) == 1);
    val[i] = macro.p[9+i];
  }
  for (int i=0; i<3; i++) xp[i] = macro.p[15+i];
  rad = macro.p[18];

  mesh(domain(typ,id)).setBoundaryConditions(fix);

  mesh(domain(typ,id)).setTmpPresDisp(flag,val,xp,rad);
*/
//--------------------------------------------------------------------------------------------------
  return 0;  
}

