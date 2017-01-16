#include "Macro.h"
#include "DomainTree.h"
#include "IsogeometricFEM.h"
//#include "Plot.h"


extern DomainTree domain;
//extern Plot       plot;


int macro224(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "bfgs";
    macro.type = "chen";
    macro.what = " BFGS solver ";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;

    macro.db.selectDomain();


    macro.db.frameButtonBox();

    macro.db.addToggleButton("line-srch/no-line-srch",true);

    macro.db.addTextField("Niter = ",10);
    
    macro.db.addTextField("STol = ",0.8,4);
    
    macro.db.addTextField("Tol = ",0.0001,6);

    macro.db.frameRadioBox();


    // and other stuff

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int  type, id, Niter;
  double  stol, tol;
  bool flag1;

  type  = roundToInt(macro.p[0]);
  id    = roundToInt(macro.p[1]) - 1;
  flag1 = (roundToInt(macro.p[2]) == 1); // line search flag
  Niter = roundToInt(macro.p[3]);
  stol   = macro.p[4];
  tol   = macro.p[5];
  
  //cout << '\t' << flag1 << '\t' << Niter << '\t' << stol << '\t' << tol << endl;

  isogeometricFEM(domain(type,id)).BFGSsolver(Niter, flag1, stol, tol);


//--------------------------------------------------------------------------------------------------
  return 0;  
}

