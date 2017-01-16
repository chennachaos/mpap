
#include "Macro.h"
#include "DomainTree.h"


extern DomainTree domain;



int macro98(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "flex";
    macro.type = "wulf";
    macro.what = "lifting line and twist analysis of a flexible wing";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;

    macro.db.selectDomain();
 
    macro.db.addRadioBox("*analyse system for (a0,q)",
                          "analyse system for (CL,q)",
                          "optimise system for (CL,q,tol)",
                          "analyse D(v) for (rho,W,v0,v1,n)",
                          "optimise D(v) for (rho,W,v0,v1,n,tol)");
 
    macro.db.addRadioBox("*coupled","decoupled");

    macro.db.addTextField("",0.,10);
    macro.db.addTextField("",0.,10);
    macro.db.addTextField("",0.,10);
    macro.db.addTextField("",0.,10);
    macro.db.addTextField("",0.,10);
    macro.db.addTextField("",0.,10);

    macro.db.stringTextField("output file name: ","out.dat",20);

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int type, id, task, i;
  
  double par[6];

  bool flag;

  type =  roundToInt(macro.p[0]);
  id   =  roundToInt(macro.p[1]) - 1;
  task =  roundToInt(macro.p[2]);
  flag = (roundToInt(macro.p[3]) == 1);
  for (i=0; i<6; i++) par[i] = macro.p[4+i];

  domain(type,id).doForFlexibleWing(task,flag,par,macro.strg);

//--------------------------------------------------------------------------------------------------
  return 0;  
}

