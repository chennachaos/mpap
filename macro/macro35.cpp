
#include "Macro.h"
#include "Definitions.h"
#include "DomainTree.h"


extern DomainTree domain;


using namespace std;


int macro35(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "dtim";
    macro.type = "outp";
    macro.what = "print domain time statistics";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;

    macro.db.selectDomain();

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int typ, id;

  typ = roundToInt(macro.p[0]);
  id  = roundToInt(macro.p[1]) - 1;

  domain(typ,id).printComputerTime();

  cout << "\n";

//--------------------------------------------------------------------------------------------------
  return 0;  
}

