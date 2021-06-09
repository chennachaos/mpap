
#include "Macro.h"


int macro6(Macro &macro)
{
  if (!macro)
  {
    macro.name = "ndif";
    macro.type = "ctrl";
    macro.what = "close if statement";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;

    return 0;
  }
//--------------------------------------------------------------------------------------------------

  // nothing to be done here! 


//--------------------------------------------------------------------------------------------------
  return 0;
}

