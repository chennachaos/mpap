
#include "Macro.h"
#include "MacroQueue.h"


extern MacroQueue macroQueue;


int macro39(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "stop";
    macro.type = "ctrl";
    macro.what = "stop macro queue execution";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;
    macro.sensitivity[PRE]   = true;

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  return macroQueue.macCmd.n + 1;

//--------------------------------------------------------------------------------------------------
}

