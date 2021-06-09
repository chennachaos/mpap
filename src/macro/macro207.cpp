#include "Macro.h"
#include "DomainTree.h"
#include "HBSplineBase.h"
#include "HBSplineFEM.h"
#include "HBSplineCutFEM.h"


extern DomainTree domain;


int macro207(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "fsi";
    macro.type = "chen";
    macro.what = " FSI schemes ";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;

    macro.db.selectDomain();

    macro.db.frameButtonBox();

    macro.db.addRadioBox("*StagForcePred","StagDispPred","MonoFixedPointFP","MonoFixedPointDP");

    macro.db.frameButtonBox();

    macro.db.addTextField("Niter = ",10);
    
    macro.db.addTextField("Tol = ",0.0001,6);

    macro.db.frameRadioBox();

    // and other stuff

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int  domType, id, Niter, soltype;
  double  tol;
  bool flag1;

  domType  = roundToInt(macro.p[0]);
  id       = roundToInt(macro.p[1]) - 1;
  soltype  = roundToInt(macro.p[2]);
  Niter    = roundToInt(macro.p[3]);
  tol      = macro.p[4];

  //cout << " domType = " << domType << endl;
  //cout << " id = " << id << endl;
  //cout << " soltype = " << soltype << endl;

  if(domType == 27)
  {
    if(soltype == 1)
      hbsplineFEM(domain(domType,id)).fsi_staggered_force_predictor(Niter, tol);
    else if(soltype == 2)
      hbsplineFEM(domain(domType,id)).fsi_staggered_displacement_predictor(Niter, tol);
    else if(soltype == 3)
      hbsplineFEM(domain(domType,id)).fsi_monolithic_fixedpoint_forcePred(Niter, tol);
    else if(soltype == 4)
      hbsplineFEM(domain(domType,id)).fsi_monolithic_fixedpoint_dispPred(Niter, tol);
    else
    {
      cout << " macro230 ... FSI schemes ... undefined for this domain " << endl;
    }
  }
  else if(domType == 28)
  {
    if(soltype == 1)
      hbsplineCutFEM(domain(domType,id)).fsi_staggered_force_predictor(Niter, tol);
    else if(soltype == 2)
      hbsplineCutFEM(domain(domType,id)).fsi_staggered_displacement_predictor(Niter, tol);
    else if(soltype == 3)
      hbsplineCutFEM(domain(domType,id)).fsi_monolithic_fixedpoint_forcePred(Niter, tol);
    else if(soltype == 4)
      hbsplineCutFEM(domain(domType,id)).fsi_monolithic_fixedpoint_dispPred(Niter, tol);
    else
    {
      cout << " macro230 ... FSI schemes ... undefined for this domain " << endl;
    }
  }

  //--------------------------------------------------------------------------------------------------

  return 0;  
}

