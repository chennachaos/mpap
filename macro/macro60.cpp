
#include "Macro.h"
#include "DomainTree.h"
//#include "Aeroplane.h"


extern DomainTree domain;


int macro60(Macro &macro)
{
/*
  if (!macro) 
  { 
    macro.name = "slup";
    macro.type = "wulf";
    macro.what = "solve and update the flight simulator";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;
    
    macro.db.selectDomain();

    //macro.db.addTextField("      element ",1);
    //macro.db.addTextField("  perturbation",0.0001,8);
    //macro.db.addTextField("digits [total]",12);
    //macro.db.addTextField("digits [ < 1 ]",5);

    //macro.db.addRadioBox("g-format","*f-format");

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int type, id;
  
  type = roundToInt(macro.p[0]);
  id   = roundToInt(macro.p[1]) - 1;

  if (!isAeroplane(domain(type,id))) return 0; 

  domain(type,id).factoriseSolveAndUpdate();
*/
//--------------------------------------------------------------------------------------------------
  return 0;  
}

