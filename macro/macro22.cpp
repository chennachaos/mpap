
#include "Macro.h"
#include "Definitions.h"
//#include "Plot.h"
#include "DomainTree.h"
//#include "Mesh.h"


extern DomainTree domain;
////extern Plot       plot;


using namespace std;


int macro22(Macro &macro)
{
/*
  if (!macro) 
  { 
    macro.name = "node";
    macro.type = "plot";
    macro.what = "plot nodes";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;
    macro.sensitivity[PRE]   = true;
    
    macro.db.selectDomain();

    //macro.db.stringList("outl","*surf","mesh");
    
    macro.db.addList(COLOURS_WHITE);

    macro.db.addToggleButton("numbers", true);

    macro.db.frameButtonBox();
	    
    macro.db.addRadioBox("initial","*current");
    
    macro.db.addRadioBox("*all","masters","dependents");
    
    macro.db.frameRadioBox();

    return 0;    
  }
//--------------------------------------------------------------------------------------------------

  int    type, id, col, num, which;
  double xmn[3], xmx[3];
  bool   def;
  
  type  = roundToInt(macro.p[0]);
  id    = roundToInt(macro.p[1]) - 1;
  col   = roundToInt(macro.p[2]);
  num   = roundToInt(macro.p[3]);
  def   = (roundToInt(macro.p[4]) == 2);
  which = roundToInt(macro.p[5]);
  
  plot.setColour(col-1);
  
  if (!plot) 
  {
     domain(type,id).findMinMaxX(xmn,xmx,def);

     plot.fit(xmn,xmx,domain(type,id).ndm);
  }
  
  domain(type,id).plotNodes(num,which,def);
*/
//--------------------------------------------------------------------------------------------------
  return 0;  
}

