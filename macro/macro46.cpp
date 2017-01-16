
#include "Macro.h"
#include "Definitions.h"
#include "Plot.h"
#include "DomainTree.h"
#include "Mesh.h"


extern DomainTree domain;
extern Plot       plot;


using namespace std;


int macro46(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "mupd";
    macro.type = "anly";
    macro.what = "ALE mesh update";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;
    macro.sensitivity[PRE]   = true;
   
    macro.db.selectDomain();

    macro.db.addRadioBox("*standard increment cutting","strict increment cutting");
    
    return 0;    
  }
//--------------------------------------------------------------------------------------------------

  int    type, id, cutType;
  
  type    = roundToInt(macro.p[0]);
  id      = roundToInt(macro.p[1]) - 1;
  cutType = roundToInt(macro.p[2]);
 
  domain(type,id).updateMesh(cutType);

//--------------------------------------------------------------------------------------------------
  return 0;  
}

