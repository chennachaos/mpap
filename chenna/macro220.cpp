
#include "Macro.h"
#include "DomainTree.h"
#include "PlotVTK.h"

extern PlotVTK  plotvtk;


int macro220(Macro &macro)
{

  if (!macro) 
  { 

    macro.name = "wrps";
    macro.type = "vtk";
    macro.what = "Write post script file";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;

//    macro.db.selectDomain();

    macro.db.stringTextField("file name ","vtkpsimage",40);

    macro.db.frameButtonBox();

    macro.db.addRadioBox("*PS","EPS","PDF","TEX","SVG");

    // and other stuff

    return 0;	

  }
//--------------------------------------------------------------------------------------------------

  int  index=0;

  index  = macro.p[0] - 1;

//cout << " index " << index << endl;

  plotvtk.VtkGL2PSExporter(macro.strg, index);


//--------------------------------------------------------------------------------------------------

  return 0;  
}

