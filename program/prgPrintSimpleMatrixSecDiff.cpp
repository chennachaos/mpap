
#include <iostream>

#include "MathBasic.h"



void prgPrintSimpleMatrixSecDiff(double* c1, double* c2, int nrow, int ncol, 
                                 int row0, int col0, int dd, 
		                 int dig, int dig2, bool gfrmt, 
	                         int indent, bool numbers)
{
  int i, j, r = row0, c = col0, dr = dd, dc = dd;

  if (r + dr - 1 > nrow) r = nrow + 1 - dr;
  if (c + dc - 1 > ncol) c = ncol + 1 - dc;
  if (r < 1) r = 1;
  if (c < 1) c = 1;
  if (r + dr - 1 > nrow) dr = nrow + 1 - r;
  if (c + dc - 1 > ncol) dc = ncol + 1 - c;
  
  char frmt[30], numberFrmt[30];
  
  if (gfrmt) sprintf(frmt," %%%d.%dg",dig,dig2);
  else       sprintf(frmt," %%%d.%df",dig,dig2);

  if (numbers)
  {
    sprintf(numberFrmt," %%%dd",dig);
    for (j=0; j<indent+dig+1; j++) printf(" ");
    for (j=c; j<c+dc; j++) printf(numberFrmt,j); printf("\n");
  }

  for (i=r; i<r+dr; i++)
  { 
    for (j=0; j<indent; j++) printf(" ");
    if (numbers) printf(numberFrmt,i);
    for (j=c; j<c+dc; j++)  printf(frmt,c1[(j-1)*nrow+i-1]-c2[(j-1)*nrow+i-1]); printf("\n");
  }
  printf("\n");

  return;
}

