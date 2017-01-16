#include "Macro.h"
#include "DomainTree.h"
#include "IsogeometricFEM.h"
//#include "Plot.h"


extern DomainTree domain;
//extern Plot       plot;


int macro227(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "sol3";
    macro.type = "chen";
    macro.what = " LM/Convection solver ";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;

    macro.db.selectDomain();


    macro.db.frameButtonBox();

    macro.db.addRadioBox("*LM","Convection");

    macro.db.frameButtonBox();

    macro.db.addTextField("Niter = ",10);
    
    macro.db.addTextField("STol = ",0.8,4);
    
    macro.db.addTextField("Tol = ",0.0001,6);

    macro.db.frameRadioBox();


    // and other stuff

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int  type, id, Niter, soltype;
  double  stol, tol;
  bool flag1;

  type  = roundToInt(macro.p[0]);
  id    = roundToInt(macro.p[1]) - 1;
  soltype = roundToInt(macro.p[2]);
  Niter = roundToInt(macro.p[3]);
  stol   = macro.p[4];
  tol   = macro.p[5];
  
  //cout << '\t' << flag1 << '\t' << Niter << '\t' << stol << '\t' << tol << endl;
  
  if(soltype == 1)
    isogeometricFEM(domain(type,id)).LMSolver(Niter, stol, tol, 1.0);
  else
    isogeometricFEM(domain(type,id)).ConvectionSolver(Niter, stol, tol, 1.0);


//--------------------------------------------------------------------------------------------------
  return 0;  
}

