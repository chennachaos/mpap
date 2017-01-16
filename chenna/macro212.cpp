
#include "Macro.h"
#include "DomainTree.h"
#include "IsogeometricFEM.h"
#include "HBSplineFEM.h"
#include "Plot.h"


extern DomainTree domain;
extern MpapTime mpapTime;


int macro212(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "sol2";
    macro.type = "chen";
    macro.what = "calculate residual, solve and update";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;

    macro.db.selectDomain();

    macro.db.addRadioBox("don't show residuals","show some residuals","*show all residuals");


    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int  type, id, printRes;

  type = roundToInt(macro.p[0]);
  id   = roundToInt(macro.p[1]) - 1;
  printRes = roundToInt(macro.p[2]);


 //hbsplineFEM(domain(type,id)).solveTimeStep();

//--------------------------------------------------------------------------------------------------
  return 0;  
}

