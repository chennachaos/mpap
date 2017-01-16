
#include "Macro.h"
#include "Plot.h"
#include "FunctionsEssGrp.h"


extern Plot plot;


int macro17(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "wipe";
    macro.type = "plot";
    macro.what = "wipe the plot area";

    macro.sensitivity[INTER]  = true;
    macro.sensitivity[BATCH]  = true;
    macro.sensitivity[NOPROJ] = true;
    macro.sensitivity[PRE]    = true;

    macro.db.addToggleButton("unset view settings", false);
    
    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  if (roundToInt(macro.p[0]) == 1) plot.resetCoor();

  essGrpWipe();
  
//--------------------------------------------------------------------------------------------------
  return 0;  
}

