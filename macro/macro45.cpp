
#include "Macro.h"
#include "Definitions.h"
//#include "Plot.h"
#include "DomainTree.h"
//#include "Mesh.h"


extern DomainTree domain;
////extern Plot       plot;


using namespace std;


int macro45(Macro &macro)
{
/*
  if (!macro) 
  { 
    macro.name = "fixd";
    macro.type = "plot";
    macro.what = "plot fixed mesh d.o.f.";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;
    macro.sensitivity[PRE]   = true;
    
    macro.db.selectDomain();

    macro.db.addList(COLOURS_YELLOW);

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
  
  if (!domain(type,id).isALE() && !isMesh(domain(type,id))) 
  {
    COUT << "Error! This domain does not have any ALE capabilities.\n\n"; return 0;
  }
  
  plot.setColour(col-1);
  
  if (!plot) 
  {
     domain(type,id).findMinMaxX(xmn,xmx,def);

     plot.fit(xmn,xmx,domain(type,id).ndm);
  }
  
  domain(type,id).plotFixed(def);
*/
//--------------------------------------------------------------------------------------------------
  return 0;  
}

