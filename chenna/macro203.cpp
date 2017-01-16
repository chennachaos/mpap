#include "Macro.h"
#include "DomainTree.h"
#include "IsogeometricFEM.h"
#include "Plot.h"


extern DomainTree domain;
extern Plot       plot;


int macro203(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "sflg";
    macro.type = "chen";
    macro.what = " Set Flags for Output operations ";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;

    macro.db.selectDomain();


    macro.db.frameButtonBox();

    macro.db.addRadioBox("*RSYS","defUndefFlag");
    
    macro.db.frameRadioBox();


    macro.db.frameButtonBox();

    macro.db.addTextField("value = ",0);

    macro.db.frameRadioBox();


    // and other stuff

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int  type, id, rsys, val;

  type = roundToInt(macro.p[0]);
  id   = roundToInt(macro.p[1]) - 1;
  rsys = roundToInt(macro.p[2]);
  val  = roundToInt(macro.p[3]);

  isogeometricFEM(domain(type,id)).setDifferentFlags(rsys, val);

  
  


//--------------------------------------------------------------------------------------------------
  return 0;  
}

