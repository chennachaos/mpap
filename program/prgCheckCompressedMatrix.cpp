
#include "MathVector.h"
#include "FunctionsProgram.h"


void prgCheckCompressedMatrix(int nRow, int nCol,
                              VectorArray<double> &x,
                              VectorArray<int> &row,
                              VectorArray<int> &col,
                              VectorArray<int> &compr)
{
  char fct[] = "prgCheckCompressedMatrix";

  int i, j;

  bool diagOK;

  if (nRow != nCol) prgError(1,fct,"nRow != nCol");

  if (row.n != col.n) prgError(2,fct,"row.n != col.n");

  if (x.n != col.n) prgError(3,fct,"x.n != col.n");

  if (compr.n < nRow) prgError(4,fct,"compr.n < nRow");

  // check for invalid rows and columns  

  for (i=0; i<row.n; i++)
  {
    if (row[i] < 1 || row[i] > nRow) prgError(5,fct,"invalid row number!");
    if (col[i] < 1 || col[i] > nCol) prgError(6,fct,"invalid column number!");
  }

  // check row compression

  VectorArray<int> cnt;

  cnt.setDim(nRow);
  cnt.zero();

  for (i=0; i<nRow; i++)
  {
    diagOK = false;

    if (row[compr[i]-1] != i+1) prgError(7,fct,"row[compr[i]-1] != i+1");

    if (row[compr[i]-1] == col[compr[i]-1]) diagOK = true;

    if (compr[i+1] - compr[i] < 1) prgError(8,fct,"empty row!");

    cnt[col[compr[i]-1]-1]++;

    for (j=compr[i]; j<compr[i+1]-1; j++)
    {
      if (row[j] != i+1) prgError(9,fct,"row[j] != i+1");

      if (col[j] == row[j]) diagOK = true;

      if (col[j] <= col[j-1]) prgError(10,fct,"col[j] <= col[j-1]");

      cnt[col[j]-1]++;
    }

    if (!diagOK) prgError(11,fct,"!diagOK");
  }

  for (i=0; i<nRow; i++) if (cnt[i] < 1) prgError(12,fct,"empty column!");

  COUT << "no errors found in row compressed matrix pattern!\n\n";

  return;
}









