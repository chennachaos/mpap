
#include "Macro.h"
#include "FunctionsMaterial.h"


int macro37(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "mat";
    macro.type = "diff";
    macro.what = "verify solid material tangent tensor";

    macro.sensitivity[NOPROJ] = true;
    macro.sensitivity[PRE]    = true;
    macro.sensitivity[INTER]  = true;
    macro.sensitivity[BATCH]  = true;
    
    macro.db.addTextField("  perturbation",0.0001,8);
    macro.db.addTextField("digits [total]",12);
    macro.db.addTextField("digits [ < 1 ]",5);

    macro.db.addRadioBox("g-format","*f-format");

    macro.db.stringTextField("data file ","matDiff.dat");
    
    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int    type, id, el, gp, dig, dig2;

  bool   gfrmt;

  double ddd;
  
  ddd    = macro.p[0];
  dig    = roundToInt(macro.p[1]);
  dig2   = roundToInt(macro.p[2]);
  gfrmt  = (roundToInt(macro.p[3]) == 1);
  
  matDiffTangentTest(ddd,macro.strg.asCharArray(),dig,dig2,gfrmt);

//--------------------------------------------------------------------------------------------------
  return 0;  
}

