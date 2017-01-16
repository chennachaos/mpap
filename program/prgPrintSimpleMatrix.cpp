
#include <iostream>

#include "MathBasic.h"



void prgPrintSimpleMatrix(double* c, int nrow, int ncol, int dig, int dig2, bool gfrmt, 
		          int indent, bool numbers)
{
  int i, j;

  char frmt[30], numberFrmt[30];
  
  if (gfrmt) sprintf(frmt," %%%d.%dg",dig,dig2);
  else       sprintf(frmt," %%%d.%df",dig,dig2);

  if (numbers)
  {
    sprintf(numberFrmt," %%%dd",dig);
    for (j=0; j<indent+dig+1; j++) printf(" ");
    for (j=1; j<ncol+1; j++)    printf(numberFrmt,j); printf("\n");
  }
 
  for (i=0; i<nrow; i++)
  { 
    for (j=0; j<indent; j++) printf(" ");
    if (numbers) printf(numberFrmt,i+1);
    for (j=0; j<ncol; j++)  printf(frmt,c[j*nrow+i]); printf("\n");
  }
  printf("\n");

  return;
}

