
#include <iostream>

#include "WorkSpace.h"
#include "Debug.h"


using namespace std;


WorkSpace::WorkSpace(void)
{
  dblSpc = NULL;
  intSpc = NULL;

  dimDbl = 0;
  dimInt = 0;
	
  if (debug) cout << " WorkSpace constructor\n\n";  

  return;
}




WorkSpace::~WorkSpace()
{
  if (dblSpc != NULL) delete [] dblSpc;
  if (intSpc != NULL) delete [] intSpc;

  if (debug) cout << " WorkSpace destructor\n\n";  

  return;
}




void WorkSpace::expandDbl(int n)
{
  if (n <= dimDbl) return;

  double *tmp = new double [n];

  if (dblSpc != NULL) delete [] dblSpc;  dblSpc = tmp;

  dimDbl = n;
	
  return;
}




void WorkSpace::expandInt(int n)
{
  if (n <= dimInt) return;

  int *tmp = new int [n];

  if (intSpc != NULL) delete [] intSpc;  intSpc = tmp;

  dimInt = n;
	
  return;
}

