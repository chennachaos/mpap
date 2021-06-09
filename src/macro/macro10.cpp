
#include "Macro.h"
#include "DomainTree.h"
//#include "Mesh.h"


using namespace std;

extern DomainTree domain;


int macro10(Macro &macro)
{
  if (!macro)
  {
    macro.name = "info";
    macro.type = "anly";
    macro.what = "print basic information";

    macro.sensitivity[BATCH] = true;
    macro.sensitivity[INTER] = true;
    macro.sensitivity[PRE]   = true;

    macro.db.selectDomain();

    return 0;
  }
//--------------------------------------------------------------------------------------------------

  int type, id;

  type = roundToInt(macro.p[0]);
  id   = roundToInt(macro.p[1]) - 1;

  cout << "\n          " << domain[type].type << " (" << id+1 << ")\n          ===========";
  for (int i=0; i<domain[type].type.length(); i++) cout << "="; cout << "\n";

  domain(type,id).printInfo();

//--------------------------------------------------------------------------------------------------
  return 0;  
}

