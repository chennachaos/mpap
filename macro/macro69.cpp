
#include "Macro.h"
#include "Mesh.h"
#include "DomainTree.h"
#include "Definitions.h"


extern DomainTree domain;


int macro69(Macro &macro)
{
  char *domTypeName[] = DOMAIN_TYPE_NAMES;

  if (!macro) 
  { 
    macro.name = "infl";
    macro.type = "prep";
    macro.what = "generate input file; append if file exists";

    macro.sensitivity[PRE] = true;

    macro.db.selectDomain();

    macro.db.addList(domTypeName);

    macro.db.addToggleButton("append to existing file",1);

    macro.db.addToggleButton("load interface discretisation",0);

    macro.db.addToggleButton("time_functions and run_control",0);

    macro.db.stringTextField("input file name ","",20);

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int typ, id, domType;

  typ = roundToInt(macro.p[0]);
  id  = roundToInt(macro.p[1]) - 1;

  domType = roundToInt(macro.p[2]);

  mesh(domain(typ,id)).writeInputFile(macro.strg.asCharArray(),domType,
                                      (roundToInt(macro.p[3]) == 1),
                                      (roundToInt(macro.p[4]) == 1),
                                      (roundToInt(macro.p[5]) == 1));

//--------------------------------------------------------------------------------------------------
  return 0;  
}

