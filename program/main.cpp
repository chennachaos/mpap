//#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET

#include <iostream>

#include "Definitions.h"
#include "FunctionsProgram.h"

#include "Global.h"
#include "MyString.h"
#include "FunctionsEssGrp.h"


#include "UnixGlobal.h"

#include "petscmat.h"
#include "conditionalOStream.h"


int main(int argc, char **argv)
{
  //PetscErrorCode ierr;

  //PetscInitialize(NULL,NULL,(char *)0,NULL);

  PetscInitialize(NULL, NULL, "petsc_options.dat", NULL);

  PetscInt  this_mpi_proc;
  MPI_Comm_rank(MPI_COMM_WORLD, &this_mpi_proc);

  ConditionalOStream  p0cout(std::cout,  (this_mpi_proc == 0) );


  int i;
  MyString tmp;  
  char *var, *cmdLineArg [] = COMMAND_LINE_ARGUMENTS;  

  computerTime.go("prepare");

  // interprete command line arguments 

  for(i=1; i<argc; i++)
  { 
    tmp = argv[i];
    switch (tmp.which(cmdLineArg))
    {  
      case  0: noGUI           = true; break;
      case  1: batch           = true; 
               noGUI           = true; break;
      case  2: debug           = true; break;
      case  3: wulf            = true; break;
      case  4: sony            = true; break;
      case  5: deniz           = true; break;
      case  6: keep            = true; break;
      case  7: aspectRatioCorr = true; break;
      case  8: readIfileInfo   = true; break;
      case  9: lastProj        = true; break;
      case 10: test            = true; break;
      case 11: break;
      case 12: break;
    }
  }

  if (test) 
  {
    prgTest();

    return 0;
  }

  // load Macros  

  prgLoadMacros();

  // print MPAP2 logo

  //prgPrintLogo();

  // noGUI is made as default option in order 
  // to be able to run the code in parallel
  
  //noGUI = true;

  // start application loop

  if(noGUI) 
  {
    // no graphical user interface

    if (!prgReadFileNames()) goto quit;

    prgReadFile();
 
    prgExecNextRunCtrlCmd();
  }
  else
  {
    // graphical user interface

    essGrpStartGUI(argc, argv);
  }

  quit:

  //delete [] mpapVar;

  //std::cout << "   =====================================================\n";
  //std::cout << "    mpap2 aborted. \n\n";
  p0cout << " \n\n\n   program is successful... \n\n\n " << endl;

  //PetscErrorCode ierr;

  PetscFinalize(); //CHKERRQ(ierr);

  return 0;
}
//

