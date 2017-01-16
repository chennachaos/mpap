
#include "Macro.h"
#include "DomainTree.h"
#include "DomainTypeEnum.h"


extern DomainTree domain;


int macro101(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "sogs";
    macro.type = "mart";
    macro.what = "SOlves one Gauss-Seidel iteration";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  if (domain.nDomainOfType(INTERFACEGS) != 1) 
  {
    cout << "   this only works for exactly one InterfaceGS!\n\n";
	return 0;
  }

  int typ = INTERFACEGS, id = 0;;

  domain(typ,id).solveGS();
  
//--------------------------------------------------------------------------------------------------
  return 0;  
}

