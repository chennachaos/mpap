
#include "Macro.h"
#include "DomainTree.h"
#include "Mesh.h"


extern DomainTree domain;


int macro34(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "proc";
    macro.type = "outp";
    macro.what = "nodal data post processing";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;
    macro.sensitivity[PRE]   = true;
    
    macro.db.selectDomain();

    macro.db.addRadioBox("*u1*u1 + u2*u2", "sqrt(u1*u1 + u2*u2)");
  
    macro.db.frameRadioBox();

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int type  = roundToInt(macro.p[0]), 
      id    = roundToInt(macro.p[1]) - 1, 
      selct = roundToInt(macro.p[2]), 
      i;

  Mesh *dom = (Mesh*) &(domain(type,id));
  
  int numnp = dom->numnp, ndf = dom->ndf;
  
  double *outp = &(dom->outp[0]), 
	 *u    = &(dom->u(1,1));

  for (i=0; i<numnp; i++)
  {
    switch (selct)
    {
      case  1: outp[i] = u[0+i*ndf]*u[0+i*ndf] + u[1+i*ndf]*u[1+i*ndf]; break;
	      
      case  2: outp[i] = sqrt(u[0+i*ndf]*u[0+i*ndf] + u[1+i*ndf]*u[1+i*ndf]); break;

      default: break;
    }
  }
  
//--------------------------------------------------------------------------------------------------
  return 0;  
}

