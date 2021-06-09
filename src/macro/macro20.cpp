
#include "Macro.h"
#include "DomainTree.h"


extern DomainTree domain;


int macro20(Macro &macro)
{
  if (!macro)
  {
    macro.name = "wmsh";
    macro.type = "prep";
    macro.what = "write mesh data to file";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;
    macro.sensitivity[PRE]   = true;

    macro.db.selectDomain();

    macro.db.addRadioBox("*initial conf.","deformed conf.");

    macro.db.stringTextField("file name ","mesh",40);

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int  type, id;
  bool defm;

  type = roundToInt(macro.p[0]);
  id   = roundToInt(macro.p[1]) - 1;
  defm = (roundToInt(macro.p[2]) == 1);

  domain(type,id).writeMeshToFile(macro.strg,defm);

//--------------------------------------------------------------------------------------------------
  return 0;  
}

