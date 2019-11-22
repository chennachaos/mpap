
#include "Macro.h"
#include "DomainTree.h"
#include "HBSplineFEM.h"
#include "HBSplineCutFEM.h"


extern DomainTree domain;
//extern MpapTime mpapTime;


int macro225(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "refi";
    macro.type = "chen";
    macro.what = "perform adaptive refinement";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;

    macro.db.selectDomain();

    macro.db.addTextField(" tol = ",0.01,5);

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int  type, id; 
  double  tol;

  type = roundToInt(macro.p[0]);
  id   = roundToInt(macro.p[1]) - 1;
  tol  = macro.p[2];
  
  hbsplineFEM(domain(type,id)).performAdaptiveRefinement(tol);

//--------------------------------------------------------------------------------------------------
  return 0;  
}

