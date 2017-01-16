
#include <iostream>

#include "SolverWorkSpace.h"
#include "Debug.h"


using namespace std;


SolverWorkSpace::SolverWorkSpace(void)
{
  dbl = NULL;

  dimDbl = 0;
	
  if (debug) cout << " SolverWorkSpace constructor\n\n";  

  return;
}




SolverWorkSpace::~SolverWorkSpace()
{
  if (dbl != NULL) delete [] dbl;

  if (debug) cout << " SolverWorkSpace destructor\n\n";  

  return;
}




void SolverWorkSpace::expand(int n)
{
  if (n <= dimDbl) return;
  
  double *tmp = new double [n];

  if (dbl != NULL) delete [] dbl;  dbl = tmp;

  dimDbl = n;
	
  return;
}

