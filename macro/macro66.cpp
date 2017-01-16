
#include "Macro.h"
#include "Mesh.h"
#include "DomainTree.h"


extern DomainTree domain;


int macro66(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "sns";
    macro.type = "prep";
    macro.what = "simple node selection";

    macro.sensitivity[PRE]   = true;
    macro.sensitivity[INTER] = true;

    macro.db.selectDomain();

    macro.db.addRadioBox("*all","inside box","inside cylinder","inside sphere");

    macro.db.addToggleButton("print nodelist",1);
    macro.db.addToggleButton("deselect",0);
    macro.db.addToggleButton("on deformed mesh",0);

    macro.db.nextButtonBox();
    macro.db.addTextField("",0);
    macro.db.addTextField("",0);
    macro.db.addTextField("",0);
    macro.db.addTextField("",0);
    macro.db.addTextField("",0);
    macro.db.addTextField("",0);
    macro.db.addTextField("",0);

    macro.db.setButtonColumnDim(3,3);

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int typ, id;
 
  typ = roundToInt(macro.p[0]);
  id  = roundToInt(macro.p[1]) - 1;

  SelectNode &selectNode = mesh(domain(typ,id)).selectNode;

  selectNode.flag = roundToInt(macro.p[2])+2;
  selectNode.prnt = (roundToInt(macro.p[3]) == 1);
  selectNode.dslt = (roundToInt(macro.p[4]) == 1);
  selectNode.defm = (roundToInt(macro.p[5]) == 1);

  for (int i=0; i<7; i++) selectNode.x[i] = macro.p[6+i];

  mesh(domain(typ,id)).simpleNodeSelection();

//--------------------------------------------------------------------------------------------------
  return 0;  
}

