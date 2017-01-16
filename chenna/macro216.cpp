
#include "Macro.h"
#include "DomainTree.h"
#include "PlotVTK.h"

extern PlotVTK  plotvtk;


int macro216(Macro &macro)
{

  if (!macro) 
  { 

    macro.name = "clr";
    macro.type = "vtk";
    macro.what = "clear in  vtk window ";

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

  plotvtk.clearWindow();


//--------------------------------------------------------------------------------------------------

  return 0;  
}

