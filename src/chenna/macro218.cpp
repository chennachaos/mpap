
#include "Macro.h"
#include "DomainTree.h"
//#include "PlotVTK.h"

//extern PlotVTK  plotvtk;


int macro218(Macro &macro)
{

  if (!macro) 
  { 

    macro.name = "wrt";
    macro.type = "vtk";
    macro.what = "Write uGrid to file";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;

//    macro.db.selectDomain();

    macro.db.stringTextField("file name ","vtkoutput.vtu",40);

    macro.db.frameButtonBox();

    // and other stuff

    return 0;	

  }
//--------------------------------------------------------------------------------------------------

  int  type, id;

  type = roundToInt(macro.p[0]);
  id   = roundToInt(macro.p[1]) - 1;

  //plotvtk.write2file(macro.strg);

//--------------------------------------------------------------------------------------------------

  return 0;  
}

