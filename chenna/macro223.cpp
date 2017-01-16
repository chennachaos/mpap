
#include "Macro.h"
#include "DomainTree.h"
#include "PlotVTK.h"

extern PlotVTK  plotvtk;


int macro223(Macro &macro)
{

  if (!macro) 
  { 

    macro.name = "refl";
    macro.type = "vtk";
    macro.what = "Reflect the data about";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;

//    macro.db.frameButtonBox();

    macro.db.addRadioBox("SetPlaneToXMin", "SetPlaneToYMin", "SetPlaneToZMin", "SetPlaneToXMax", "SetPlaneToYMax", "SetPlaneToZMax","*SetPlaneToX","SetPlaneToY","SetPlaneToZ");

    // and other stuff

    return 0;	

  }
//--------------------------------------------------------------------------------------------------

  int  index=0;

  index  = macro.p[0] - 1;

  plotvtk.VtkReflect(index);


//--------------------------------------------------------------------------------------------------

  return 0;  
}

