
#include "Macro.h"
#include "Definitions.h"
//#include "Plot.h"
#include "DomainTree.h"
//#include "Mesh.h"


extern DomainTree domain;
////extern Plot       plot;


using namespace std;


int macro23(Macro &macro)
{
/*
  if (!macro) 
  { 
    macro.name = "elem";
    macro.type = "plot";
    macro.what = "plot element numbers";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;
    macro.sensitivity[PRE]   = true;
    
    macro.db.selectDomain();

    macro.db.addList(COLOURS_GREEN);

    macro.db.addRadioBox("initial conf.","*current conf.");
    
    macro.db.frameRadioBox();

    return 0;    
  }
//--------------------------------------------------------------------------------------------------

  int    type, id, col;
  double xmn[3], xmx[3];
  bool   def;
  
  type = roundToInt(macro.p[0]);
  id   = roundToInt(macro.p[1]) - 1;
  col  = roundToInt(macro.p[2]);
  def  = (roundToInt(macro.p[3]) == 2);
  
  plot.setColour(col-1);
  
  if (!plot) 
  {
     domain(type,id).findMinMaxX(xmn,xmx,def);

     plot.fit(xmn,xmx,domain(type,id).ndm);
  }
  
  domain(type,id).plotElemNum(def);
*/
//--------------------------------------------------------------------------------------------------
  return 0;  
}

