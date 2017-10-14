

#include "Macro.h"
#include "DomainTree.h"
//#include "PlotVTK.h"

//extern PlotVTK  plotvtk;


int macro222(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "axes";
    macro.type = "vtk";
    macro.what = "Add Axes to the vtk window";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;


    macro.db.frameButtonBox();

    macro.db.addRadioBox("*Add axes","Remove axes");
    
    macro.db.frameRadioBox();
    
    macro.db.addTextField(" x = ",0.2,3);
    macro.db.addTextField(" y = ",0.2,3);
    macro.db.addTextField(" z = ",0.2,3);

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  bool  flag;
  double u1, v1, w1;

  flag = (macro.p[0] == 1);

  u1 = macro.p[1];
  v1 = macro.p[2];
  w1 = macro.p[3];
  
  //if(flag)  
    //plotvtk.addAxes(u1, v1, w1);
  //else
    //plotvtk.removeAxes();

//--------------------------------------------------------------------------------------------------
  return 0;  
}

