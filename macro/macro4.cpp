
#include "Macro.h"
#include "MacroQueue.h"
#include "Loop.h"


extern MacroQueue macroQueue;


int macro4(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "xlop";
    macro.type = "ctrl";
    macro.what = "exit active loop";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;
    
    return 0;  
  }
//--------------------------------------------------------------------------------------------------

  int  i = roundToInt(macro.p[0]), maxIter;
  Loop *loop = &(macroQueue.loop[i]);

  loop->cnt = 0;

  return loop->end + 2;

//--------------------------------------------------------------------------------------------------
}

