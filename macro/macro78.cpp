
#include "Macro.h"
#include "DomainTree.h"
#include "SolverSparse.h"
#include "FiniteElementBVP.h"


extern DomainTree domain;


int macro78(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "matx";
    macro.type = "anly";
    macro.what = "plot global system matrix pattern";

    macro.sensitivity[BATCH] = true;
    macro.sensitivity[INTER] = true;

    macro.db.selectDomain();

    macro.db.stringTextField("eps file name (optional)","",40);

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int  type, id;

  char *fileName;

  type = roundToInt(macro.p[0]);
  id   = roundToInt(macro.p[1]) - 1;

  if (!!macro.strg) fileName = macro.strg.asCharArray(); else fileName = NULL;

  domain(type,id).plotSolverMatrixPattern(fileName);

//--------------------------------------------------------------------------------------------------
  return 0;  
}

