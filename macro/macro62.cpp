
#include "Macro.h"
#include "SolverTime.h"


extern SolverTime solverTime;



int macro62(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "stim";
    macro.type = "outp";
    macro.what = "print linear solver statistics";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;

    // and other stuff

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  solverTime.print();

//--------------------------------------------------------------------------------------------------
  return 0;  
}

