
#include <iostream>

#include "RunControl.h"
#include "FunctionsProgram.h"
#include "Debug.h"
#include "MyString.h"
#include "MacroCommand.h"
#include "MacroList.h"


extern RunControl runCtrl;
extern MacroList  macro;


void prgGetMacroCommand(MacroCommand &macCmd, int &macj)
{
  char *mode [] = RUN_MODE_PROMPTS;
	
  MyString macStrg;
 
  int cmdi = runCtrl.cmd[runCtrl.act];
  
  while (1)
  { 
    if (cmdi == 1) 
      { std::cout << mode[runCtrl.mode] << ":" << runCtrl.macCnt+1 << "> "; 
        macStrg.input().stripToMin(); }
    else
      { if (macj >= runCtrl.batch[-cmdi].n) 
          prgError(1,"prgGetMacroCommand","batch macros do not conclude with 'end'.");
        macStrg = runCtrl.batch[-cmdi][macj++]; }

   
    if (prgStringToMacroCmd(macCmd,macStrg)) break;
  }

  runCtrl.macCnt++;
  
  return;
}


