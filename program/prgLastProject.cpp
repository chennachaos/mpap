
#include <iostream>
#include <fstream>

#include "Files.h"
#include "FunctionsProgram.h"
#include "Definitions.h"
#include "Debug.h"
#include "MyString.h"


extern Files files;


bool prgLastProject(void)
{ 
  std::ifstream in;	
  
  // check for file "mpap.name" in current directory

  in.open(LAST_PROJECT);

  if (!in) return false;

  // read old project file names

  files.projDir.read(in);
  files.Ifile.read(in);
  files.Ofile.read(in);
  files.Tfile.read(in);
  files.Pfile.read(in);
  
  in.close();

  if (debug)
    {  std::cout << " " << files.projDir << "\n";
       std::cout << " " << files.Ifile << "\n";
       std::cout << " " << files.Ofile << "\n";
       std::cout << " " << files.Tfile << "\n";
       std::cout << " " << files.Pfile << "\n\n"; }
 

  // check for obvious errors

  if (!files.projDir.DirectoryOK()) return false;
  if (!files.Ifile.FileOK()) return false;
  if (!files.Ofile.FileOK()) return false;
  if (!files.Tfile.FileOK()) return false;
  if (!files.Pfile.FileOK()) return false;
 
  return true; 
}


