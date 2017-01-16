
#include "Macro.h"
#include "DomainTree.h"
#include "PlotVTK.h"

extern PlotVTK  plotvtk;


int macro215(Macro &macro)
{

  if (!macro) 
  { 

    macro.name = "bgcl";
    macro.type = "vtk";
    macro.what = "set background color in  vtk";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;

//    macro.db.addList(COLOURS_BLUE);

    macro.db.addTextField("R ", 0.4);

    macro.db.addTextField("G ", 0.1);

    macro.db.addTextField("B ", 0.2);

    macro.db.addTextField("A ", 0);

//    macro.db.frameButtonBox();

    // and other stuff

    return 0;	

  }
//--------------------------------------------------------------------------------------------------

  double  color[4];

  color[0] = macro.p[0];
  color[1] = macro.p[1];
  color[2] = macro.p[2];
  color[3] = macro.p[3];

  plotvtk.setBackgroundColor(color);


//--------------------------------------------------------------------------------------------------

  return 0;  
}

