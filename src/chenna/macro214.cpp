
#include "Macro.h"
#include "DomainTree.h"


int macro214(Macro &macro)
{

  if (!macro) 
  { 

    macro.name = "set";
    macro.type = "vtk";
    macro.what = "activate/deactivate vtk control";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;

//    macro.db.selectDomain();

    macro.db.frameButtonBox();

    // and other stuff

    return 0;	

  }
//--------------------------------------------------------------------------------------------------

  int  type, id;

  type = roundToInt(macro.p[0]);
  id   = roundToInt(macro.p[1]) - 1;

  //plotvtk.set();

//--------------------------------------------------------------------------------------------------

  return 0;  
}

