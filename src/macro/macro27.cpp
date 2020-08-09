
#include "Macro.h"
#include "FunctionsProgram.h"
#include "MacroQueue.h"
#include "ComputerTime.h"
#include "FunctionsEssGrp.h"
#include "RunControl.h"
#include "petscmat.h"


extern MacroQueue   macroQueue;
extern ComputerTime computerTime;
extern RunControl   runCtrl;
extern bool         noGUI;


int macro27(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "wait";
    macro.type = "ctrl";
    macro.what = "time delay or user interaction";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;
    
    macro.db.stringList("dlay","keyb","*mous");

    macro.db.addTextField(" delay [s/1000] = ",10,7);

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  MyString ch; 
  
  char *yes[] = YESS, *no[] = NOO;

  if (macro.strg == "keyb" || (macro.strg == "mous" && noGUI)) 
  {
    while (ch.which(yes)<0 && ch.which(no)<0) 
    {
      //cout << "    continue with macro execution ? (y/n) [y] ";
      PetscPrintf(MPI_COMM_WORLD, "    continue with macro execution ? (y/n) [y] ");

      ch.free().append("y").inputKeepIfReturn();
    }
    //p0cout << "\n";

    if (ch.which(yes) < 0) return macroQueue.macCmd.n + 1; else return 0;
  } 
  else if (macro.strg == "mous")
  {
    PetscPrintf(MPI_COMM_WORLD, "continue with left mouse button or\n");
    PetscPrintf(MPI_COMM_WORLD, "stop with  any other mouse button!\n\n");

    /*
    essGrpSetSensAllButDrawingArea(false);
    runCtrl.fixStatus(PRESSMOUSE);

    if (essGrpWaitForMouseButtonPressed() != 1) 
    {
      essGrpSetSensAllButDrawingArea(true);
      runCtrl.freeFixedStatus();
      return macroQueue.macCmd.n + 1;
    }
    else
    {
      essGrpSetSensAllButDrawingArea(true);
      runCtrl.freeFixedStatus();
      return 0;
    }
    */
  }
  
  computerTime.sleep(macro.p[0]);

//--------------------------------------------------------------------------------------------------
  return 0;  
}

