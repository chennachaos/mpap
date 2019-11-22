
#include "Macro.h"


void optimizeFlightPath(void);


int macro61(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "path";
    macro.type = "wulf";
    macro.what = "optimize flight path";

    // and other stuff

    macro.sensitivity[NOPROJ] = true;

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  //optimizeFlightPath();

//--------------------------------------------------------------------------------------------------
  return 0;  
}

