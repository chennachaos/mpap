
#include "Macro.h"
#include "DomainTree.h"
#include "Plot.h"
#include "BeamSection.h"


extern DomainTree domain;
extern Plot plot;


int macro55(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "sect";
    macro.type = "anly";
    macro.what = "analyse shear flow in 3D beam";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;
    macro.sensitivity[PRE]   = true;
    
    macro.db.selectDomain();

    macro.db.addList(COLOURS_BLUE);                   //  2
                                                         
    macro.db.addToggleButton("section outline",1);    //  3
    macro.db.addToggleButton("geom. centre CA",0);    //  4
    macro.db.addToggleButton("shear centre CS",0);    //  5
    macro.db.addToggleButton("I,J & prcpl.ax.",0);    //  6
    macro.db.addToggleButton("plot S in A"    ,0);    //  7
    macro.db.addToggleButton("plot M in CA"   ,0);    //  8
    macro.db.addToggleButton("booms >> skin area",0); //  9
    macro.db.addToggleButton("plot distribut.",0);    // 10

    macro.db.frameButtonBox();
    macro.db.frameRadioBox();

    macro.db.addRadioBox("*q (Sy' in CS)","q (Sz' in CS)","q (T')","q (S in A & T)",  // 11
                         "tau (S in A & T)","sig (M and N)","sig* (N,S,M,T)");

    macro.db.addToggleButton("with labels",0); // 12
    macro.db.addToggleButton("wipe first",0);  // 13
    macro.db.addToggleButton("zoom out",0);    // 14
    macro.db.addToggleButton("reset view",0);  // 15
    macro.db.addTextField(" scale ",1.);       // 16
    macro.db.addTextField("interv.",16);       // 17

    macro.db.setButtonColumnDim(20,20);

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int    type, id, col, distr, n;
  double xmn[2], xmx[2], scl, fact = 1.2;
  
  type = roundToInt(macro.p[0]);
  id   = roundToInt(macro.p[1]) - 1;
  col  = roundToInt(macro.p[2]);

  if (!isBeamSection(domain(type,id))) 
    { cout << "          'sect' is for domain type BEAMSECTION only!\n\n"; return 0; }

  bool section = (roundToInt(macro.p[3]) == 1),
       CA      = (roundToInt(macro.p[4]) == 1),
       CS      = (roundToInt(macro.p[5]) == 1),
       geom    = (roundToInt(macro.p[6]) == 1),
       SinA    = (roundToInt(macro.p[7]) == 1),
       MinCA   = (roundToInt(macro.p[8]) == 1), 
       bOnly   = (roundToInt(macro.p[9]) == 1), labels;

  distr = roundToInt(macro.p[11]);

  if (roundToInt(macro.p[10]) != 1) distr = 0;
  labels  = (roundToInt(macro.p[12]) == 1);

  scl = macro.p[16];
 
  n = roundToInt(macro.p[17]);

  if (roundToInt(macro.p[13]) == 1) plot.wipe();
  if (roundToInt(macro.p[14]) == 1) 
  {
    plot.x0Des[0] -= (fact - 1.) * .5 * plot.dDes[0];
    plot.x0Des[1] -= (fact - 1.) * .5 * plot.dDes[1];
    plot.dDes[0]  *= fact;
    plot.dDes[1]  *= fact;
    plot.adjustToNewSize();
  }

  plot.setColour(col-1);
  
  if (!plot || roundToInt(macro.p[15]) == 1) 
  {
     domain(type,id).findMinMaxX(xmn,xmx);

     plot.fit(xmn,xmx,domain(type,id).ndm);
  }

  domain(type,id).doForSection(section,CA,CS,bOnly,geom,SinA,MinCA,labels,distr,scl,n);

//--------------------------------------------------------------------------------------------------
  return 0;  
}

