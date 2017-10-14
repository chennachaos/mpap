
#include "FunctionsMaterial.h"
#include "FunctionsProgram.h"
#include "Definitions.h"


int matdim_(int *matId)
{
  int materialDim[] = MATERIAL_DIMENSIONS;

  //std::cout << *matId << " -> " << materialDim[*matId-1] <<"\n";
  
  if (*matId < 0) prgError(1,"matDim","invalid matId!");

  int dim = materialDim[*matId-1];

  //cout << *matId << "->" << dim << "\n";

  if (dim < 1 || dim > 3) 
    prgError(2,"matDim","invalid matId or MATERIAL_DIMENSIONS not up to date!");
  
  return dim;
}


