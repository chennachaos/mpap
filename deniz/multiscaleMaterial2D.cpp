
#include "DomainTree.h"
#include "Element.h"
#include "FunctionsProgram.h"
#include "FunctionsDeniz.h"


extern DomainTree domain;


void multiscalematerial2d_(double *matData, double *F, double *stre, double *cc, 
		           int *finite, int *gp, int *error, Element *elm)

{
  int id;
  
  if(elm!=NULL) id = elm->GpDatId[*gp-1];

  else  prgError(7,"multiscalematerial2d","multiscalematerial2d called without elm!!!");

  *error=domain(MICROCELL,id).microSolve2D(matData, F, stre, cc, (*finite==1));

  return;
}

