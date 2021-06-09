
#include "Macro.h"
#include "DomainTree.h"


extern DomainTree domain;



int macro13(Macro &macro)
{
  if (!macro)
  {
    macro.name = "solv";
    macro.type = "anly";
    macro.what = "select and initialise solver";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;

    macro.db.selectDomain();

    macro.db.addList("*MA41","PARDISO(sym)","PARDISO(unsym)");

    macro.db.addTextField("parameter 1: ",1,5);
    macro.db.addTextField("parameter 2: ",1,5);
    macro.db.addTextField("parameter 3: ",1,5);

    macro.db.addToggleButton("checkIO",false);

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int    type, id, slv, parm[10];

  bool   cIO;

  type    = roundToInt(macro.p[0]);
  id      = roundToInt(macro.p[1]) - 1;
  slv     = roundToInt(macro.p[2]);
  parm[0] = roundToInt(macro.p[3]);
  parm[1] = roundToInt(macro.p[4]);
  parm[2] = roundToInt(macro.p[5]);
  cIO     = (roundToInt(macro.p[6]) == 1);

 // cout << type << '\t' << id << endl;

  domain(type,id).setSolver(slv,parm,cIO);

//--------------------------------------------------------------------------------------------------
  return 0;  
}

