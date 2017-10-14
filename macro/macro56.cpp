
#include "Macro.h"
#include "DomainTree.h"
//#include "Plot.h"
//#include "BeamBending.h"


extern DomainTree domain;
////extern Plot plot;


int macro56(Macro &macro)
{
/*
  if (!macro) 
  { 
    macro.name = "bend";
    macro.type = "anly";
    macro.what = "analyse bending of 3D beam";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;
    macro.sensitivity[PRE]   = true;
    
    macro.db.selectDomain();                   // 0,1

    macro.db.addList(COLOURS_CYAN);            // 2

    macro.db.addToggleButton("system",1);      // 3
    macro.db.addToggleButton("deflection",0);  // 4
    macro.db.addToggleButton("rotation",0);    // 5
    macro.db.addToggleButton("moment",0);      // 6
    macro.db.addToggleButton("shear force",0); // 7
    macro.db.addToggleButton("XZ-plane",1);    // 8
    macro.db.addToggleButton("XY-plane",1);    // 9

    macro.db.frameButtonBox();
    macro.db.nextButtonBox();

    macro.db.addToggleButton("with labels",0); // 10
    macro.db.addToggleButton("wipe first",0);  // 11
    macro.db.addToggleButton("reset view",0);  // 12
    macro.db.addTextField("interv.",16);       // 13

    macro.db.setButtonColumnDim(20,20);

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int    type, id, col, distr, n;
  double xmn[2], xmx[2];
  
  type = roundToInt(macro.p[0]);
  id   = roundToInt(macro.p[1]) - 1;
  col  = roundToInt(macro.p[2]);

  if (!isBeamBending(domain(type,id))) 
    { cout << "          'bend' is for domain type BEAMBENDING only!\n\n"; return 0; }

  bool syst  = (roundToInt(macro.p[ 3]) == 1),
       defl  = (roundToInt(macro.p[ 4]) == 1),
       rota  = (roundToInt(macro.p[ 5]) == 1),
       bmom  = (roundToInt(macro.p[ 6]) == 1),
       shfo  = (roundToInt(macro.p[ 7]) == 1),
       XZ    = (roundToInt(macro.p[ 8]) == 1),
       XY    = (roundToInt(macro.p[ 9]) == 1),
       label = (roundToInt(macro.p[10]) == 1);

  if (!XZ && !XY)
  {
    COUT << "choose XY, XZ or both!\n\n";
    return 0;
  }

  n = roundToInt(macro.p[13]);

  if (roundToInt(macro.p[11]) == 1) plot.wipe();

  plot.setColour(col-1);
  
  if (!plot || roundToInt(macro.p[12]) == 1)
  {
     domain(type,id).findMinMaxXBB(syst,defl,rota,bmom,shfo,XZ,XY,xmn,xmx);

     plot.fit(xmn,xmx,domain(type,id).ndm,0.);
  }

  domain(type,id).doForBending(syst,defl,rota,bmom,shfo,XZ,XY,label,n);
*/
//--------------------------------------------------------------------------------------------------
  return 0;  
}

