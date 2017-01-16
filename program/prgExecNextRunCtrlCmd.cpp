

#include "RunControl.h"
#include "FunctionsProgram.h"
#include "MacroCommand.h"
#include "MacroList.h"
#include "RunModeEnum.h"
#include "MyStatusEnum.h"
#include "MacroQueue.h"



extern MacroQueue macroQueue;
extern RunControl runCtrl;
extern MacroList  macro;
extern bool       noGUI;



void prgExecNextRunCtrlCmd(void)
{
  if (++runCtrl.act > runCtrl.cmd.n - 1) return;

  MacroCommand macCmd;

  int macj, cmdi = runCtrl.cmd[runCtrl.act];

  switch (cmdi)
  {
     default: // BATCH and INTER : runCtrl.cmd[i] < 2
    	     
              if (cmdi == 1) { runCtrl.newMode(INTER); runCtrl.newStatus(INTERACTIVE); }
	      else           { runCtrl.newMode(BATCH); runCtrl.newStatus(BATCHMODE);   }
	    
              if (cmdi == 1 && !noGUI) return;  // return control to GUI 
		
              macj = 0;
	
              prgGetMacroCommand(macCmd,macj);
	
              while (macro[macCmd.ii].name != "end")
              {
                //cout << macro[macCmd.ii].name << "\n";

                macroQueue.append(macCmd); 

                prgGetMacroCommand(macCmd,macj);  
              }

              macroQueue.append(macCmd);
  
              return;

     //case  2: ........ break;

     //case  3: ........ break;
  }

  std::cout << " Should this really have happened ?\n\n";
 
  prgExecNextRunCtrlCmd();
 
  return;
}

