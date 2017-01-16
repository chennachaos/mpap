
#include "Macro.h"
#include "Definitions.h"
#include "Plot.h"
#include "DomainTree.h"
#include "Mesh.h"


extern DomainTree domain;
extern Plot       plot;


using namespace std;


int macro26(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "reac";
    macro.type = "plot";
    macro.what = "plot reaction forces";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;
    macro.sensitivity[PRE]   = true;
    
    macro.db.selectDomain();

    macro.db.addList(COLOURS_LIGHTBLUE);

    macro.db.addRadioBox("initial conf.","*current conf.");
    
    macro.db.frameRadioBox();
    
    macro.db.addTextField(" maximum force = ",-1,6);

    macro.db.addLabel(" set < 0 for autom. scaling ");
    
    macro.db.frameButtonBox();
    
    return 0;    
  }
//--------------------------------------------------------------------------------------------------

  int    type, id, col;
  double xmn[3], xmx[3], scl;
  bool   def;
  
  type = roundToInt(macro.p[0]);
  id   = roundToInt(macro.p[1]) - 1;
  col  = roundToInt(macro.p[2]);
  def  = (roundToInt(macro.p[3]) == 2);
  scl  = macro.p[4];
  
  plot.setColour(col-1);
  
  if (!plot) 
  {
     domain(type,id).findMinMaxX(xmn,xmx,def);

     plot.fit(xmn,xmx,domain(type,id).ndm);
  }
  
  domain(type,id).plotReac(scl,def);

//--------------------------------------------------------------------------------------------------
  return 0;  
}

