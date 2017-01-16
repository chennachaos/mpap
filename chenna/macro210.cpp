
#include "Macro.h"
#include "DomainTree.h"
#include "IsogeometricFEM.h"
#include "Plot.h"
#include "HBSplineFEM.h"

extern DomainTree domain;
extern Plot       plot;

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
  rr = roundToInt(macro.p[2]);
  geom = macro.p[3];
  index = macro.p[4];
  patchnum = macro.p[5];

  hbsplineFEM(domain(type,id)).writeReadResult(rr, macro.strg);
  
  


//--------------------------------------------------------------------------------------------------
  return 0;  
}

