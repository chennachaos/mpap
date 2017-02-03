
#include <iostream>

#include "FunctionsProgram.h"
#include "Files.h"
#include "Definitions.h"
#include "Debug.h"
#include "MyString.h"
#include "petscmat.h"


extern Files files;
extern bool  lastProj;

using namespace std;


bool prgReadFileNames()
{
  MyString answ;
  char *no [] = NOO;
  char *yes [] = YESS;
  char *quit [] = QUIT;

  bool ok = false;
  
  
  if(prgLastProject())
  {  
    PetscPrintf(MPI_COMM_WORLD, "       project directory   : %s   \n", files.projDir.asCharArray());
    PetscPrintf(MPI_COMM_WORLD, "       input file name     : %s   \n", files.Ifile.asCharArray());
    PetscPrintf(MPI_COMM_WORLD, "       output file name    : %s   \n", files.Ofile.asCharArray());
    PetscPrintf(MPI_COMM_WORLD, "       time plot file name : %s   \n", files.Tfile.asCharArray());
    PetscPrintf(MPI_COMM_WORLD, "       eps plot file name  : %s \n\n", files.Pfile.asCharArray());

    if (!lastProj)
    {
      do
      {
        PetscPrintf(MPI_COMM_WORLD, "       Are file names correct ? (y/n/x) "); answ.input();
      } 
      while ( answ.which(no)<0 && answ.which(yes)<0 && answ.which(quit)<0 );
      PetscPrintf(MPI_COMM_WORLD, "\n\n");

      if(answ.which(quit) >= 0)
        return false;

      if(answ.which(yes)  >= 0 && prgFileExist(files.projDir,files.Ifile))
        ok = true;
      else if(answ.which(no) < 0) 
        PetscPrintf(MPI_COMM_WORLD, "      No such project directory or input file!\n\n");
    }
    else
    {
      if(prgFileExist(files.projDir,files.Ifile))
        return true;
      else
      {
        PetscPrintf(MPI_COMM_WORLD, "       No such project directory or input file!\n");
        return false;
      }
    }
  }

  
  while(!ok)
  {
    PetscPrintf(MPI_COMM_WORLD, "       project directory   : "); files.projDir.input();
    do
    {
      PetscPrintf(MPI_COMM_WORLD, "       input file name     : "); files.Ifile.input();
    }
    while (files.Ifile.length()<1);

    files.Ofile.free().append(files.Ifile)[0] = 'O';
    files.Tfile.free().append(files.Ifile)[0] = 'T';
    files.Pfile.free().append(files.Ifile)[0] = 'P';

    PetscPrintf(MPI_COMM_WORLD, "       output file name    : %s ", files.Ofile.asCharArray()); files.Ofile.inputKeepIfReturn();
    PetscPrintf(MPI_COMM_WORLD, "       time plot file name : %s ", files.Tfile.asCharArray()); files.Tfile.inputKeepIfReturn();
    PetscPrintf(MPI_COMM_WORLD, "       eps plot file name  : %s ", files.Pfile.asCharArray()); files.Pfile.inputKeepIfReturn(); 
    PetscPrintf(MPI_COMM_WORLD, "\n");

    do
    {
      PetscPrintf(MPI_COMM_WORLD, "       Are file names correct ? (y/n/x) "); answ.input();
    }
    while ( answ.which(no)<0 && answ.which(yes)<0 && answ.which(quit)<0 );
    PetscPrintf(MPI_COMM_WORLD, "\n\n");

    if(answ.which(quit) >= 0)
      return false;

    if(answ.which(yes)  >= 0 && prgFileExist(files.projDir,files.Ifile))
      ok = true;
    else if (answ.which(no) < 0)
      PetscPrintf(MPI_COMM_WORLD, "       No such project directory or input file!\n\n");
  }

  if(debug)
  {
    PetscPrintf(MPI_COMM_WORLD, "  %s \n ", files.projDir.asCharArray());
    PetscPrintf(MPI_COMM_WORLD, "  %s \n ", files.Ifile.asCharArray());
    PetscPrintf(MPI_COMM_WORLD, "  %s \n ", files.Ofile.asCharArray());
    PetscPrintf(MPI_COMM_WORLD, "  %s \n ", files.Tfile.asCharArray());
    PetscPrintf(MPI_COMM_WORLD, "  %s \n ", files.Pfile.asCharArray());
  }

  prgWriteLastProject();

  return true;
}


