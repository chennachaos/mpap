
#include "Macro.h"

#include "mpi.h"


int macro41(Macro &macro)
{
  int  this_mpi_proc;
  MPI_Comm_rank(MPI_COMM_WORLD, &this_mpi_proc);

  if (!macro) 
  { 
    macro.name = "prnt";
    macro.type = "outp";
    macro.what = "print string in console window";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;
    macro.sensitivity[PRE]   = true;

    macro.db.stringTextField("string to be printed: ");

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  if(this_mpi_proc == 0)
  {
    COUT << macro.strg << "\n";
    if(macro.strg.length() > 0)
      COUT << "\n";
  }

//--------------------------------------------------------------------------------------------------
  return 0;  
}

