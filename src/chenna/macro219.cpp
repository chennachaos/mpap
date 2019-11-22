
#include "Macro.h"
#include "DomainTree.h"
//#include "PlotVTK.h"

//extern PlotVTK  plotvtk;


int macro219(Macro &macro)
{

  if (!macro) 
  { 

    macro.name = "wrim";
    macro.type = "vtk";
    macro.what = "Write image file";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;
    macro.sensitivity[PRE]   = true;

//    macro.db.selectDomain();

    macro.db.stringTextField("file name","vtkoutput",40);

    macro.db.frameButtonBox();

    macro.db.addRadioBox("*PNG","JPG","TIFF","BMP","PNM","PS");

    return 0;	

  }
//--------------------------------------------------------------------------------------------------

  int  index=0;

  index  = macro.p[0] - 1;

  //plotvtk.writeImages(macro.strg, index);

//--------------------------------------------------------------------------------------------------

  return 0;  
}

