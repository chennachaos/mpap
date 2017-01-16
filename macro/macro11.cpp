
#include "Macro.h"
#include "DomainTree.h"
#include "MpapTime.h"


extern DomainTree domain;
extern MpapTime mpapTime;


int macro11(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "tang";
    macro.type = "anly";
    macro.what = "calculate stiffness, residual, solve and update";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;

    macro.db.selectDomain();

    macro.db.addRadioBox("don't show residuals","show some residuals","*show all residuals");

    return 0;    
  }
//--------------------------------------------------------------------------------------------------

  int type     = roundToInt(macro.p[0]),
      id       = roundToInt(macro.p[1]) - 1,
      printRes = roundToInt(macro.p[2]);

  if (!domain(type,id).solverOK) 
    {  COUT << "use 'solv' to initialise a solver first!\n\n";  return 0;  }
  
  if (!mpapTime.dtOK) 
    {  COUT << "use 'dt' to set the time step size first!\n\n"; return 0;  }

  if (domain(type,id).calcStiffnessAndResidual(printRes) != 0) return 0;

  if (domain(type,id).converged()) return 0;

  if (domain(type,id).factoriseSolveAndUpdate()  != 0) return 0;

  domain(type,id).updateIterStep();
  
//--------------------------------------------------------------------------------------------------
  return 0;  
}

