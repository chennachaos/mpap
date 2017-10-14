
#include "Macro.h"
#include "DomainTree.h"


extern DomainTree domain;


int macro47(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "ale";
    macro.type = "diff";
    macro.what = "verify ALE element stiffness";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;
    
    macro.db.selectDomain();

    macro.db.addTextField("      element ",1);
    macro.db.addTextField("  perturbation",0.0001,8);
    macro.db.addTextField("digits [total]",12);
    macro.db.addTextField("digits [ < 1 ]",5);

    macro.db.addRadioBox("g-format","*f-format");
    
    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int    type, id, el, dig, dig2;

  bool   gfrmt;

  double ddd;
  
  type   = roundToInt(macro.p[0]);
  id     = roundToInt(macro.p[1]) - 1;
  el     = roundToInt(macro.p[2]);
  ddd    = macro.p[3];
  dig    = roundToInt(macro.p[4]);
  dig2   = roundToInt(macro.p[5]);
  gfrmt  = (roundToInt(macro.p[6]) == 1);
  
  domain(type,id).elementDiffStiffTestMesh(ddd,el,dig,dig2,gfrmt);

//--------------------------------------------------------------------------------------------------
  return 0;  
}

