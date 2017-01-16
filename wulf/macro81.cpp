
#include "Macro.h"
#include "DomainTree.h"
#include "MicroCellWulf.h"
#include "FunctionsProgram.h"


extern DomainTree domain;


int macro81(Macro &macro)
{
  if (!macro) 
  { 
    /*macro.name = "defm";
    macro.type = "wulf";
    macro.what = "get uDep data from strain (2D 'MicroCell' only)";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;
    
    macro.db.selectDomain();

    macro.db.addTextField("eps_11 = ",0.);
    macro.db.addTextField("eps_22 = ",0.);
    macro.db.addTextField("eps_12 = ",0.);
   
    macro.db.frameButtonBox();
*/
    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int    type, id;
  double strain[10];
  
  type      = roundToInt(macro.p[0]);
  id        = roundToInt(macro.p[1]) - 1;
  strain[0] = macro.p[2];
  strain[1] = macro.p[3];
  strain[2] = macro.p[4];
 
  if (!isMicroCellWulf(domain(type,id))) 
    prgWarning(1,"macro51","this works only for domains of type MicroCellWulf!");
  
  domain(type,id).strainToBoundaryDisplacement(strain);

//--------------------------------------------------------------------------------------------------
  return 0;  
}

