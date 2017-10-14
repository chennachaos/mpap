
#include "Macro.h"
#include "MpapTime.h"
#include "MacroQueue.h"
#include "Loop.h"


extern MpapTime   mpapTime;
extern MacroQueue macroQueue;


int macro2(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "loop";
    macro.type = "ctrl";
    macro.what = "begin macro loop";

    macro.db.stringTextField("loop name :","",30);
     
    macro.db.addTextField("iterations or time :", 10, 8);
     
    macro.db.addRadioBox("*count iterations","stop at max time");
     
    macro.sensitivity[INTER] = true; 
    macro.sensitivity[BATCH] = true; 
    macro.sensitivity[PRE]   = true; 
     
    return 0;    
  }
//--------------------------------------------------------------------------------------------------

  int    i = roundToInt(macro.p[0]), maxIter;
  Loop   *loop = &(macroQueue.loop[i]);
  double maxTime;
  
  if (roundToInt(macro.p[2]) == 2)
  {
    maxTime = macro.p[1];
    if (maxTime <= mpapTime.cur) return loop->end + 2;
  }
  else
  { 
    maxIter = roundToInt(macro.p[1]);
    if (loop->cnt > maxIter) { loop->cnt = 1; return loop->end + 2; }
  }

//--------------------------------------------------------------------------------------------------
  return 0;  
}

