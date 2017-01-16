
#include "Macro.h"
#include "Plot.h"
#include "FunctionsEssGrp.h"


extern Plot plot;


int macro31(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "draw";
    macro.type = "plot";
    macro.what = "refresh drawing area";

    macro.sensitivity[BATCH] = true;

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

   essGrpCopyPixmap(); 

//--------------------------------------------------------------------------------------------------
  return 0;  
}

