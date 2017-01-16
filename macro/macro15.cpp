
#include "Macro.h"
#include "MpapTime.h"
#include "DomainTree.h"
#include "FunctionsProgram.h"

#include "petscmat.h"


extern DomainTree domain;
extern MpapTime mpapTime;


int macro15(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "dt";
    macro.type = "anly";
    macro.what = "define time step size";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;

    macro.db.stringList("*init","fact","adap","cut","stck");
    
    macro.db.addTextField("dt  / fact / iter / = ",0.1,6);
    macro.db.addTextField("max(dt)//opt.iter / = ",0,6);
    
    return 0;    
  }
//--------------------------------------------------------------------------------------------------

  if (macro.strg == "init")
  {  
    mpapTime.dtOK     = true;
    mpapTime.dt       = fabs(macro.p[0]);
    mpapTime.dtMax    = fabs(macro.p[1]);
    mpapTime.stack.free();
    mpapTime.stack.append(mpapTime.dt);
    if (mpapTime.dtMax < 1.e-15) mpapTime.dtMax = mpapTime.dt;
  }
  else if (!mpapTime.dtOK) 
  {
    prgWarning(1,"macro15","dt macro ignored! initialise dt first!"); return 0;
  }	  
  else if (macro.strg == "fact")
  {
    mpapTime.dt *= fabs(macro.p[0]);
  }
  else if (macro.strg == "adap")
  {
    mpapTime.dt *= pow(0.5,(macro.p[0] - macro.p[1]));
  }
  else if (macro.strg == "cut")
  {
    mpapTime.cut();
  }
  else if (macro.strg == "stck")
  {
    mpapTime.stck();
  }
  if(mpapTime.dt > mpapTime.dtMax)
    mpapTime.dt = mpapTime.dtMax;
  
  for(int i=0; i<domain.ndom; i++)
    domain(i).setTimeParam(); 

  PetscPrintf(MPI_COMM_WORLD, "  time increment dt = %10.6f \n\n", mpapTime.dt);

//--------------------------------------------------------------------------------------------------
  return 0;  
}

