
#include "Macro.h"
#include "ComputerTime.h"
#include "FunctionsProgram.h"


extern ComputerTime computerTime;



int macro79(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "ctim";
    macro.type = "outp";
    macro.what = "computer time stop watch";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;

    macro.db.stringTextField("identifier name","STOPTIME",40);

    macro.db.addToggleButton("write to TFile",false);

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  double dt;

  char tmp[30], *strg;

  if (macro.strg == "\0") macro.strg = "STOPTIME";

  strg = macro.strg.asCharArray();

  if (!computerTime.go(strg,false))
  {
    dt = computerTime.stop(strg);

    COUT << "computer time " << strg << " = " << dt << " sec\n\n";

    if (roundToInt(macro.p[0]) == 1)
    {
      sprintf(tmp,"%12.5g",dt);

      prgWriteToTFile(tmp);
    }
  }

//--------------------------------------------------------------------------------------------------
  return 0;  
}

