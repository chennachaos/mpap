

#include "Macro.h"
#include "DomainTree.h"
#include "HBSplineFEM.h"
#include "HBSplineCutFEM.h"
#include "StandardFEM.h"

extern DomainTree domain;


int macro205(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "res1";
    macro.type = "chen";
    macro.what = "print disps & stresses etc";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;
    macro.sensitivity[PRE]   = true;
    
    macro.db.selectDomain();

    macro.db.addTextField(" patch = ",1);
    macro.db.addTextField(" u = ",0,7);
    macro.db.addTextField(" v = ",0,7);
    macro.db.addTextField(" w = ",0,7);

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int type, id, patchnum;
  double u1, v1, w1;
  
  type = roundToInt(macro.p[0]);
  id   = roundToInt(macro.p[1]) - 1;

  patchnum = roundToInt(macro.p[2]) - 1;
  u1 = macro.p[3];
  v1 = macro.p[4];
  w1 = macro.p[5];
  
  //if(type == 26)
    //isogeometricFEM(domain(type,id)).printResultAtPoint(patchnum, u1, v1, w1);
 
  if(type == 27)
    hbsplineFEM(domain(type,id)).printResultAtPoint(patchnum, u1, v1, w1);

//--------------------------------------------------------------------------------------------------
  return 0;  
}

