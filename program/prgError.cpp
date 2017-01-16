
#include <iostream>

#include "FunctionsProgram.h"


void prgError(int id, char* fname, char* msg)
{
  std::cout << " ERROR (" << id << ") in " << fname << ": " << msg << "\n\n";

  exit(0);
}




void prgerror_(int *id, char* fname, char* msg)
{
  std::cout << " ERROR (" << *id << ") in " << fname << ": " << msg << "\n\n";

  exit(0);
}





