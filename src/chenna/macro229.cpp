#include "Macro.h"
#include "DomainTree.h"
#include "HBSplineFEM.h"
#include "HBSplineCutFEM.h"
<<<<<<< HEAD:chenna/macro229.cpp
//#include "Plot.h"


extern DomainTree domain;
//extern Plot       plot;
=======


extern DomainTree domain;
>>>>>>> collabchandan:src/chenna/macro229.cpp


int macro229(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "lfdr";
    macro.type = "chen";
    macro.what = " Compute total force on immersed body ";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;

    macro.db.selectDomain();

    //macro.db.addRadioBox("*alg1","alg2");

    macro.db.addTextField("IB id ", 1);

    macro.db.frameRadioBox();

    //macro.db.nextButtonBox();

    //macro.db.addTextField("u0 ", 0.1, 5);
    //macro.db.addTextField("u1 ", 0.3, 5);

    //macro.db.frameButtonBox();

    //macro.db.nextButtonBox();

    //macro.db.addTextField("v0 ", 0.3, 5);
    //macro.db.addTextField("v1 ", 0.7, 5);

    //macro.db.frameButtonBox();



    // and other stuff

    return 0;   
  }
//--------------------------------------------------------------------------------------------------

  int  type, id, index;

  type  = roundToInt(macro.p[0]);
  id    = roundToInt(macro.p[1]) - 1;

  index = roundToInt(macro.p[2]);
  index -= 1;

  //cout << index << '\t' << nPs << '\t' << u0 << '\t' << u1 << '\t' << v0 << '\t' << v1 << endl;
  //if(type == 26)
    //cout << " isogeometricFEM(domain(type,id)).computeLiftAndDrag()  .... not implemented yet " << endl;
 
  if(type == 27)
    hbsplineFEM(domain(type,id)).computeTotalForce(index);

  if(type == 28)
    hbsplineCutFEM(domain(type,id)).computeTotalForce(index);

//--------------------------------------------------------------------------------------------------
  return 0;  
}

