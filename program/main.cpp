
#include <iostream>
#include "Definitions.h"
#include "FunctionsProgram.h"
#include "Global.h"
#include "MyString.h"
#include "FunctionsEssGrp.h"
#include "petscmat.h"
#include "Files.h"


extern Files files;



int main(int argc, char **argv)
{
    PetscInitialize(NULL, NULL, "petsc_options.dat", NULL);

    if(argc < 3)
    {
      // not enough input arguments. So, exit the program.

      PetscPrintf(MPI_COMM_WORLD, "  ================================== \n\n");
      PetscPrintf(MPI_COMM_WORLD, "  \n Not enough input arguments \n");
      PetscPrintf(MPI_COMM_WORLD, "  \n Program aborted... \n\n");
      PetscPrintf(MPI_COMM_WORLD, "  ================================== \n\n");

      PetscFinalize(); //CHKERRQ(ierr);

      return 1;
    }

    // load Macros  
    prgLoadMacros();

    // print MPAP2 logo
    // prgPrintLogo();
    // not doing here for parallel code. Maybe in the future, when I have more time to do "pretty" trivial stuff.

    // noGUI is made as default option in order 
    // to be able to run the code in parallel
    noGUI = true;

    // start application loop

    // set directory path and file name

    files.projDir = argv[1];
    files.Ifile   = argv[2];

    files.Ofile.free().append(files.Ifile)[0] = 'O';
    files.Tfile.free().append(files.Ifile)[0] = 'T';
    files.Pfile.free().append(files.Ifile)[0] = 'P';

    PetscPrintf(MPI_COMM_WORLD, "       project directory   : %s   \n", files.projDir.asCharArray());
    PetscPrintf(MPI_COMM_WORLD, "       input file name     : %s   \n", files.Ifile.asCharArray());
    PetscPrintf(MPI_COMM_WORLD, "       output file name    : %s   \n", files.Ofile.asCharArray());
    PetscPrintf(MPI_COMM_WORLD, "       time plot file name : %s   \n", files.Tfile.asCharArray());
    PetscPrintf(MPI_COMM_WORLD, "       eps plot file name  : %s \n\n", files.Pfile.asCharArray());

    if( prgFileExist(files.projDir,files.Ifile) )
    {
      //PetscPrintf(MPI_COMM_WORLD, "  ================================== \n\n");
      //PetscPrintf(MPI_COMM_WORLD, "  Directory and file exist \n\n");
      //PetscPrintf(MPI_COMM_WORLD, "  ================================== \n\n");

      // continue to execute the program by reading the input file
      prgReadFile();
 
      prgExecNextRunCtrlCmd();

      PetscPrintf(MPI_COMM_WORLD, "  ================================== \n\n");
      PetscPrintf(MPI_COMM_WORLD, "  Program is successful ... \n\n ");
      PetscPrintf(MPI_COMM_WORLD, "  ================================== \n\n");
    }
    else
    {
      PetscPrintf(MPI_COMM_WORLD, "  ================================== \n\n");
      PetscPrintf(MPI_COMM_WORLD, "  No such project directory or input file! \n\n");
      PetscPrintf(MPI_COMM_WORLD, "  Program has been aborted ... \n\n");
      PetscPrintf(MPI_COMM_WORLD, "  ================================== \n\n");

      PetscFinalize(); //CHKERRQ(ierr);

      return 1;
    }

    PetscFinalize(); //CHKERRQ(ierr);

    return 0;
}


