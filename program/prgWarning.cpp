
#include <iostream>

#include "FunctionsProgram.h"


void prgWarning(int id, char* fname, char* msg)
{
  std::cout << " WARNING (" << id << ") in " << fname << ": " << msg << "\n\n";

  return;
}





void prgwarning_(int *id, char* fname, char* msg)
{
  std::cout << " WARNING (" << *id << ") in " << fname << ": " << msg << "\n\n";

  return;
}





