
#include "Macro.h"
#include "DomainTree.h"


extern DomainTree domain;


int macro38(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "stif";
    macro.type = "diff";
    macro.what = "verify global stiffness matrix";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;
    
    macro.db.selectDomain();

    macro.db.addTextField("  perturbation",0.0001,8);
    macro.db.addTextField("digits [total]",11);
    macro.db.addTextField("digits [ < 1 ]",5);

    macro.db.addRadioBox("g-format","*f-format");

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int    type, id, el, gp, dig, dig2;

  bool   gfrmt;

  double ddd;
  
  type   = roundToInt(macro.p[0]);
  id     = roundToInt(macro.p[1]) - 1;
  ddd    = macro.p[2];
  dig    = roundToInt(macro.p[3]);
  dig2   = roundToInt(macro.p[4]);
  gfrmt  = (roundToInt(macro.p[5]) == 1);

  domain(type,id).globalDiffStiffTest(ddd,dig,dig2,gfrmt);

//--------------------------------------------------------------------------------------------------
  return 0;  
}

