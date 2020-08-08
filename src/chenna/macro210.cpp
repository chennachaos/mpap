
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


    macro.db.addTextField("index ", 1);

    macro.db.addTextField("stride ", 10);

    macro.db.frameButtonBox();

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int type   = roundToInt(macro.p[0]);
  int id     = roundToInt(macro.p[1]) - 1;
  int index  = roundToInt(macro.p[2]);
  int stride = roundToInt(macro.p[3]);

  //cout << id << '\t' << index << '\t' << stride << '\t' << macro.strg <<  endl;

  string  filename(macro.strg.asCharArray());

  hbsplineCutFEM(domain(type,id)).writeReadResult(index, filename, stride);

//--------------------------------------------------------------------------------------------------
  return 0;  
}

