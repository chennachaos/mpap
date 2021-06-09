
#include "Macro.h"
#include "MpapTime.h"
#include "FunctionsEssGrp.h"
#include "FunctionsProgram.h"
#include "List.h"
#include "TimeFunction.h"

#include "petscmat.h"


extern MpapTime           mpapTime;
extern List<TimeFunction> timeFunction;
extern double             globalMaxIncrement;


int macro11(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "time";
    macro.type = "anly";
    macro.what = "increment time";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;

    macro.db.stringList("*forw","back");

    return 0;
  }
//--------------------------------------------------------------------------------------------------

  if (!mpapTime.dtOK)
  {
    COUT << "    use 'dt' to set the time step size first!\n\n";
    return 0;
  }

  globalMaxIncrement = 0.;
 
  if (macro.strg == "forw")
  {
    //mpapTime.dtPrev  = mpapTime.dt;

    mpapTime.prev2 = mpapTime.prev;

    mpapTime.prev  = mpapTime.cur;

    mpapTime.cur  += mpapTime.dt;
  }
  else if (macro.strg == "back")
  {
    if (mpapTime.prev > mpapTime.prev2 + 1.e-12)
    {
      mpapTime.cur  = mpapTime.prev;
      mpapTime.prev = mpapTime.prev2;
    }
    else
    {
      prgWarning(1,"macro13","'time,back' ignored (I can only do one step back!)");
      return 0;
    }
  }
  //essGrpWriteTime();

  for(int i=0; i<timeFunction.n; i++)
    timeFunction[i].update();

  PetscPrintf(MPI_COMM_WORLD, "    time update, current time = %10.6f  \n\n", mpapTime.cur);

//--------------------------------------------------------------------------------------------------
  return 0;  
}

