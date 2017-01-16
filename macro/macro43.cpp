
#include "Macro.h"
#include "Plot.h"


extern Plot plot;


int macro43(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "post";
    macro.type = "outp";
    macro.what = "open / close postscript file";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;
    macro.sensitivity[PRE]   = true;

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  if (!plot.psOpen) 
  {
    if (!plot.psOpenFile()) return 0;

    plot.psOpen = true;
  }
  else
  {
    plot.psCloseFile();
    
    plot.psOpen = false;
  }

//--------------------------------------------------------------------------------------------------
  return 0;  
}

