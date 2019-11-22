
#include "Macro.h"
#include "Definitions.h"
#include "RunControl.h"
#include "MacroQueue.h"
#include "FunctionsEssGrp.h"


extern RunControl runCtrl;
extern MacroQueue macroQueue;


int macro64(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "test";
    macro.type = "wulf";
    macro.what = "test various stuff";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;
    macro.sensitivity[PRE]   = true;

    macro.db.selectDomain();

    macro.db.addRadioBox("test nodal data transfer for 2D meshes",
                         "select node with mouse",
                         "*mouse interaction");

    return 0;    
  }
//--------------------------------------------------------------------------------------------------

  int i, isw = roundToInt(macro.p[2]),
      x, y;
  bool inDrawingArea;

  if (isw == 1)// || isw == 2)
  {
    if (runCtrl.mode == BATCH) 
      { COUT << "mouse input not admissible in 'batch' mode!\n\n"; return 0; }

    //essGrpSetSensAllButDrawingArea(false);
   
    switch (isw)
    {
      case 1: runCtrl.fixStatus(EXECTESTMACRO); break;
      //case 2: runCtrl.fixStatus(SELECTPOINT);   break;
    }

    for (i=0; i<macro.p.n; i++) macroQueue.p[i] = macro.p[i]; 
    macroQueue.p.trunc(i);
  }
  else if (isw == 3)
  {
    //i = essGrpWasMouseButtonPressed(&x,&y,&inDrawingArea);

    //if (i != 0)  cout << "   " << i << " -> " << x << "," << y << "; " << inDrawingArea << "\n";

  }
//--------------------------------------------------------------------------------------------------
  return 0;  
}

