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
    macro.name = "erro";
    macro.type = "chen";
    macro.what = " Compute Error Norm ";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;

    macro.db.selectDomain();


    macro.db.frameButtonBox();

    macro.db.addTextField("index = ",1);
    
    macro.db.addTextField("tol = ",0.001,6);

    macro.db.frameRadioBox();


    // and other stuff

    return 0;   
  }
//--------------------------------------------------------------------------------------------------

  int  type, id, index;
  double  tol;

  type  = roundToInt(macro.p[0]);
  id    = roundToInt(macro.p[1]) - 1;

  index = roundToInt(macro.p[2]);
  tol   = macro.p[3];
  
  switch(type)
  {
    //case 26:
        //isogeometricFEM(domain(type,id)).computeElementErrors(index);
      //break;
 
    case 27:
        hbsplineFEM(domain(type,id)).computeElementErrors(index);
      break;

    case 28:
        hbsplineCutFEM(domain(type,id)).computeElementErrors(index);
      break;

    case 29:
        standardFEM(domain(type,id)).computeElementErrors(index);
      break;

    //case 30:
      //  hbscutFEMElasticity(domain(type,id)).computeElementErrors(index);
      //break;

    default :
        cout << " macro228 ... computeElementErrors() ... undefined for this domain " << endl;
      break;
  }
//--------------------------------------------------------------------------------------------------
  return 0;  
}

