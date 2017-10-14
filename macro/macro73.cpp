
#include "Macro.h"
//#include "Mesh.h"
#include "DomainTree.h"


extern DomainTree domain;


int macro73(Macro &macro)
{
/*
  if (!macro) 
  { 
    macro.name = "free";
    macro.type = "prep";
    macro.what = "set free free nodes";

    macro.sensitivity[PRE] = true;

    macro.db.selectDomain();

    macro.db.addRadioBox("*set current node selection to free nodes",
                         "unset currently selected free nodes",
                         "save current interface discretisation");
    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int typ, id, flag;

  typ = roundToInt(macro.p[0]);
  id  = roundToInt(macro.p[1]) - 1;

  flag = roundToInt(macro.p[2]);

  mesh(domain(typ,id)).setFreeNodes(flag);
*/
//--------------------------------------------------------------------------------------------------
  return 0;  
}

