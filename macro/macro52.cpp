
#include "Macro.h"
#include "DomainTree.h"
#include "Mesh.h"


extern DomainTree domain;



int macro52(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "gmsh";
    macro.type = "anly";
    macro.what = "generate mesh (call 'elsz' first!)";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;
    macro.sensitivity[PRE]   = true;

    macro.db.selectDomain();

    macro.db.addToggleButton("show mode");

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int  type, id;

  bool showFlg;

  type    = roundToInt(macro.p[0]);
  id      = roundToInt(macro.p[1]) - 1;
  showFlg = (roundToInt(macro.p[2]) == 1); 

//  if (domain(type,id).inheritsFrom(MESH)) 
//    { COUT << "this is for Meshes and derived domains only! abort gmsh.\n\n"; return 0; }

  Mesh &dom = *((Mesh*)(&domain(type,id)));

  if (dom.elemGrpToBeMeshed.n == 0)
    { COUT << "use 'elsz' to set parameters for remeshing first!\n\n"; return 0; }

  dom.remeshElemGroups(showFlg);
  
//--------------------------------------------------------------------------------------------------
  return 0;  
}

