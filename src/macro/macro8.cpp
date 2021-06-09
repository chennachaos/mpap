
#include "Macro.h"
#include "FunctionsProgram.h"
#include "RunControl.h"


extern RunControl runCtrl;


int macro8(Macro &macro)
{
  if (!macro)
  { 
    macro.name = "end";
    macro.type = "ctrl";
    macro.what = "close interactive/batch mode";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;
    macro.sensitivity[PRE]   = true;

    macro.db.addToggleButton("exit mpap2",false);

    return 0;
  }
//--------------------------------------------------------------------------------------------------

  std::cout << "          " << macro << "\n\n";

  std::cout << "          This macro will never be executed!                             \n";
  std::cout << "                                                                         \n";
  std::cout << "          The necessary actions are implemented in 'MacroQueue::append'. \n\n";

//--------------------------------------------------------------------------------------------------
  return 0;  
}

