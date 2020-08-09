
#include "Macro.h"
#include "RunControl.h"
#include "MathBasic.h"


extern RunControl runCtrl;


int macro1(Macro &macro)
{
  if (!macro) 
  { 
/*    macro.name = "mode";
    macro.type = "ctrl";
    macro.what = "set mpap2 mode";

    macro.sensitivity[INTER] = true; 
    macro.sensitivity[PRE]   = true; 
    
    macro.db.addRadioBox("pre-processor","*analysis","post-processor");
    
    macro.db.frameRadioBox();
*/
    return 0;
  }
//--------------------------------------------------------------------------------------------------
/*
  RunMode mode;
  
  if (runCtrl.mode == BATCH) return 0;
  
  if (roundToInt(macro.p[0]) <= 1) mode = PRE;
  if (roundToInt(macro.p[0]) == 2) mode = INTER;

  runCtrl.newMode(mode);

  if (runCtrl.mode == PRE)   runCtrl.newStatus(PREPRO);
  if (runCtrl.mode == INTER) runCtrl.newStatus(INTERACTIVE);
 
  macro.db.dflt[0] = macro.p[0];
*/
//--------------------------------------------------------------------------------------------------
  return 0;  
}

