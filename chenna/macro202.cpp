
#include "Macro.h"
#include "DomainTree.h"
#include "IsogeometricFEM.h"
#include "Plot.h"

extern DomainTree domain;

int macro202(Macro &macro)
{

  if (!macro) 
  { 

    macro.name = "sol";
    macro.type = "chen";
    macro.what = "solv eqns";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;

    macro.db.selectDomain();

    macro.db.frameButtonBox();

    // and other stuff

    return 0;	

  }
//--------------------------------------------------------------------------------------------------

  int  type, id;

  type = roundToInt(macro.p[0]);
  id   = roundToInt(macro.p[1]) - 1;

  //  isogeometricFEM(domain(type,id)).calcStiffnessAndResidual(2);
  //  isogeometricFEM(domain(type,id)).factoriseSolveAndUpdate();
  
    isogeometricFEM(domain(type,id)).finalsolve();


//--------------------------------------------------------------------------------------------------

  return 0;  
}

