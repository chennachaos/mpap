
#include "Macro.h"
#include "MacroQueue.h"
#include "Loop.h"


extern MacroQueue macroQueue;


int macro2(Macro &macro)
{
  if (!macro)
  {
    macro.name = "next";
    macro.type = "ctrl";
    macro.what = "close macro loop";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;
    macro.sensitivity[PRE]   = true;

    return 0;
  }
//--------------------------------------------------------------------------------------------------

  int    i = roundToInt(macro.p[0]);
  Loop   *loop = &(macroQueue.loop[i]);

  loop->cnt++;

  return loop->beg + 1;

//--------------------------------------------------------------------------------------------------
}

