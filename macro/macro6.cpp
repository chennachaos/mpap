
#include "Macro.h"
#include "MacroQueue.h"
#include "If.h"


extern MacroQueue macroQueue;


int macro6(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "else";
    macro.type = "ctrl";
    macro.what = "to be used as part of if statement";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int  i = roundToInt(macro.p[0]);
  If   *iff = &(macroQueue.iff[i]);
 
  return iff->end + 1;

//--------------------------------------------------------------------------------------------------
  return 0;  
}

