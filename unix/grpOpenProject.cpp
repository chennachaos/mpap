
#include "FunctionsProgram.h"
#include "FunctionsEssGrp.h"


void grpOpenProject(void)
{
  essGrpWriteProject();
 
  prgReadFile();

  prgExecNextRunCtrlCmd();

  return;  
}

 

