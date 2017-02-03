
#include "Macro.h"
#include "DomainTree.h"

extern DomainTree domain;


int macro211(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "mshp";
    macro.type = "chen";
    macro.what = "plot mode shape";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;

    macro.db.selectDomain();

    macro.db.frameButtonBox();

    macro.db.addTextField("Mode = ",1);

    macro.db.frameRadioBox();


    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int  type, id, printRes=1, ind;

  type = roundToInt(macro.p[0]);
  id   = roundToInt(macro.p[1]) - 1;
  ind  = roundToInt(macro.p[2]) - 1;

  //isogeometricFEM(domain(type,id)).plotModeShape(ind, 0);

//--------------------------------------------------------------------------------------------------
  return 0;  
}

