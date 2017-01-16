
#include "Macro.h"
#include "DomainTree.h"
#include "IsogeometricFEM.h"
#include "Plot.h"


extern DomainTree domain;
extern Plot       plot;


int macro206(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "pres";
    macro.type = "chen";
    macro.what = "print results";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;

    macro.db.selectDomain();


    macro.db.frameButtonBox();

    macro.db.addRadioBox("*disp","strain","stress");
    
    macro.db.frameRadioBox();

    macro.db.addTextField(" patch = ",1);

    macro.db.addTextField("DIR = ",1);

    macro.db.addTextField("param = ",0,5);

    macro.db.addTextField("incr = ",0.05,5);

    macro.db.frameButtonBox();

    // and other stuff

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  //std::cout << "          " << macro << "\n\n";

  int  type, id, dir, val1, patchnum;
  double param, incr;

  type = roundToInt(macro.p[0]);
  id   = roundToInt(macro.p[1]) - 1;
  val1 = roundToInt(macro.p[2]) - 1;
  patchnum = roundToInt(macro.p[3]) - 1;
  dir = roundToInt(macro.p[4]);
  param = macro.p[5];
  incr = macro.p[6];

  if(val1 == 0)
    isogeometricFEM(domain(type,id)).printDispsAtParameter(patchnum, dir, param, incr);
  else if(val1 == 1)
    isogeometricFEM(domain(type,id)).printStrainsAtParameter(patchnum, dir, param, incr);
  else
    isogeometricFEM(domain(type,id)).printStressesAtParameter(patchnum, dir, param, incr);


  
  


//--------------------------------------------------------------------------------------------------
  return 0;  
}

