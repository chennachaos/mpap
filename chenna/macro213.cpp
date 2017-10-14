
#include "Macro.h"
#include "DomainTree.h"


extern DomainTree domain;


int macro213(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "ccsf";
    macro.type = "chen";
    macro.what = "check continuity of 2 surfaces";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;

    macro.db.selectDomain();


    macro.db.frameButtonBox();

    macro.db.addRadioBox("*init","final");
    
    macro.db.frameRadioBox();

    macro.db.frameButtonBox();

    macro.db.addTextField(" patch1 = ",1);

    macro.db.addTextField("patch2 = ",2);

    macro.db.addTextField("dir = ",1);

    macro.db.frameButtonBox();

    // and other stuff

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int  type, id, dir, patch1, patch2, geom;

  type = roundToInt(macro.p[0]);
  id   = roundToInt(macro.p[1]) - 1;

  geom = roundToInt(macro.p[2]);
  patch1 = roundToInt(macro.p[3]) - 1;
  patch2 = roundToInt(macro.p[4]) - 1;
  dir = roundToInt(macro.p[5]);

  //isogeometricFEM(domain(type,id)).checkContSurfs(geom, patch1, patch2, dir);

//--------------------------------------------------------------------------------------------------
  return 0;  
}

