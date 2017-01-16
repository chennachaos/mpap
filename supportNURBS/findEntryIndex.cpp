#include <iostream>
#include <math.h>
#include "NurbsUtilities.h"

using namespace std;


int findEntryIndex(MatrixSparse<double>& mtx1, int& r1, int& c1)
{
  for(int i=0;i<mtx1.x.n;i++)
  {
    if(mtx1.row[i] == r1)
    {
      if(mtx1.col[i] == c1)
        return i+1;
    }
  }

  return -1;
}
