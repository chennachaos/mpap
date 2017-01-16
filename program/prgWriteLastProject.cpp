

#include <iostream>
#include <fstream>

#include "FunctionsProgram.h"
#include "Files.h"
#include "Definitions.h"
#include "MyString.h"


extern Files files;


void prgWriteLastProject(void)
{
  std::ofstream out;
	
  // write mpap2.name

  out.open(LAST_PROJECT);

  if (!out) prgError(1,"prgWriteLastProject","could not open file for writing!");
  
  out << files.projDir << "\n" 
      << files.Ifile   << "\n"
      << files.Ofile   << "\n"
      << files.Tfile   << "\n"
      << files.Pfile   << "\n";

  out.close();

  return;                                              
}


