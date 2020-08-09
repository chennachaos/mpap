
#include "Macro.h"
#include "Definitions.h"
//#include "Plot.h"
#include "DomainTree.h"
//#include "Mesh.h"
//#include "IsogeometricFEM.h"


extern DomainTree domain;
////extern Plot       plot;


using namespace std;


int macro14(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "mesh";
    macro.type = "plot";
    macro.what = "plot mesh";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;
    macro.sensitivity[PRE]   = true;
    
    macro.db.selectDomain();

    //macro.db.stringList("outl","*surf","mesh");

    macro.db.addList(COLOURS_BLUE);

    macro.db.addRadioBox("initial conf.","*current conf.");
    
    macro.db.frameRadioBox();

    macro.db.addToggleButton("backface swap (3D)",0);
    macro.db.addToggleButton("mesh outline (2D)",0);

    macro.db.frameButtonBox();
 
    return 0;    
  }
//--------------------------------------------------------------------------------------------------

  int    type, id, col;
  double xmn[3], xmx[3];
  bool   defm, flag;
  
  type = roundToInt(macro.p[0]);
  id   = roundToInt(macro.p[1]) - 1;
  col  = roundToInt(macro.p[2]);
  defm = (roundToInt(macro.p[3]) == 2);

  flag = ((roundToInt(macro.p[4]) == 1 && domain(type,id).ndm == 3) ||
          (roundToInt(macro.p[5]) == 1 && domain(type,id).ndm == 2));

  Domain &dom = domain(type,id);

  //plot.setColour(col-1);
  
  //if (!plot) 
  //{
     //dom.findMinMaxX(xmn,xmx,defm);

    // plot.fit(xmn,xmx,dom.ndm);
  //}
  
  //dom.plotMesh(defm,flag);

//--------------------------------------------------------------------------------------------------
  return 0;  
}

