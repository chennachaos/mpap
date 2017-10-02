
#include "FunctionsProgram.h"



void prgerror_(int *id, char* fname, char* msg)
{
  std::cout << " ERROR (" << *id << ") in " << fname << ": " << msg << "\n\n";
  exit(0);
}

void prgwarning_(int *id, char* fname, char* msg)
{
  std::cout << " WARNING (" << *id << ") in " << fname << ": " << msg << "\n\n";

  return;
}


