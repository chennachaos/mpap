
#include "DomainTree.h"
#include "Element.h"
#include "FunctionsDeniz.h"
#include "FunctionsProgram.h"


extern DomainTree domain;


void multiscalematerial3d_(double *matData, double *F, double *stre, double *cc, 
		           int *finite, int *gp, int *error, Element *elm)
{
  int id;
  
  if(elm!=NULL) id = elm->GpDatId[*gp-1];
  
  else  prgError(7,"multiscalematerial3d","multiscalematerial3d called without elm!!!");

  *error=domain(MICROCELL,id).microSolve3D(matData, F, stre, cc, (*finite==1));

  return;
}

