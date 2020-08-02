
#include "Macro.h"
#include "DomainTree.h"
#include "HBSplineFEM.h"
#include "HBSplineCutFEM.h"


extern DomainTree domain;

int macro210(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "wrhb";
    macro.type = "chen";
    macro.what = "write/read result";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;
    macro.sensitivity[PRE]   = true;

    macro.db.selectDomain();

    macro.db.frameButtonBox();

    macro.db.addRadioBox("*write","read");

    macro.db.frameRadioBox();

    macro.db.stringTextField("file name ","hb-result",40);

    macro.db.frameButtonBox();

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int  type, id, patchnum, geom, index, rr;

  type = roundToInt(macro.p[0]);
  id   = roundToInt(macro.p[1]) - 1;
  index = macro.p[2];

  cout << type << '\t' << id << '\t' << index << endl;

  hbsplineCutFEM(domain(type,id)).writeReadResult(index, macro.strg);

//--------------------------------------------------------------------------------------------------
  return 0;  
}

