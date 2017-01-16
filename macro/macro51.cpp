
#include "Macro.h"
#include "DomainTree.h"
#include "Plot.h"



extern DomainTree domain;
extern Plot       plot;


static const unsigned int pow2[8] = {1,2,4,8,16,32,64,128};



int macro51(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "geom";
    macro.type = "plot";
    macro.what = "plot geometry";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;
    macro.sensitivity[PRE]   = true;

    macro.db.selectDomain();

    macro.db.addList(COLOURS_RED);

    macro.db.addToggleButton("points",1);
    macro.db.addToggleButton("splines",1);
    macro.db.addToggleButton("surfaces",1);
    macro.db.addToggleButton("numbers",0);
    macro.db.addToggleButton("numbers",0);
    macro.db.addToggleButton("numbers",0);

    macro.db.setButtonColumnDim(3);

    macro.db.frameButtonBox();

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int  type, id, col;
  double xmn[3], xmx[3];
  unsigned int bitFlag = 0;

  type = roundToInt(macro.p[0]);
  id   = roundToInt(macro.p[1]) - 1;
  col  = roundToInt(macro.p[2]);

  if (roundToInt(macro.p[3])) bitFlag = bitFlag | pow2[0];
  if (roundToInt(macro.p[4])) bitFlag = bitFlag | pow2[2];
  if (roundToInt(macro.p[5])) bitFlag = bitFlag | pow2[4];
  if (roundToInt(macro.p[6])) bitFlag = bitFlag | pow2[1];
  if (roundToInt(macro.p[7])) bitFlag = bitFlag | pow2[3];
  if (roundToInt(macro.p[8])) bitFlag = bitFlag | pow2[5];
    
  plot.setColour(col-1);
  
  if (!plot) 
  {
     domain(type,id).findGeomMinMaxX(xmn,xmx);

     plot.fit(xmn,xmx,domain(type,id).ndm);
  }
  
  domain(type,id).plotGeometry(bitFlag);
 
//--------------------------------------------------------------------------------------------------
  return 0;  
}

