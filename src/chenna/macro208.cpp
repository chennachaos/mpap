
#include "Macro.h"
#include "DomainTree.h"
#include "HBSplineCutFEM.h"
#include "StandardFEM.h"
<<<<<<< HEAD:chenna/macro208.cpp
//#include "Plot.h"


extern DomainTree domain;
//extern Plot       plot;
=======


extern DomainTree domain;

>>>>>>> collabchandan:src/chenna/macro208.cpp

int macro208(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "wr2f";
    macro.type = "chen";
    macro.what = "write/read geom/mesh to a file";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;
    macro.sensitivity[PRE]   = true;

    macro.db.selectDomain();

    macro.db.frameButtonBox();

    macro.db.addRadioBox("*write","read");

    macro.db.frameRadioBox();


    macro.db.stringTextField("file name ","surfout",40);

    macro.db.frameButtonBox();

    macro.db.addRadioBox("*initialgeom","finalgeom","result");
    
//    macro.db.frameRadioBox();

    macro.db.frameButtonBox();

    macro.db.addRadioBox("*elements","control points.");
    
    macro.db.frameRadioBox();

    macro.db.addTextField("patch ",1);

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

  //if(rr == 1)
    //isogeometricFEM(domain(type,id)).writeGeomToFile(macro.strg, geom, index, patchnum);
  //if(rr == 2)
    //isogeometricFEM(domain(type,id)).readSurfaceFromFile(macro.strg, geom, patchnum);

//--------------------------------------------------------------------------------------------------
  return 0;  
}

