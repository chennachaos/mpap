
#ifndef incl_ContainerBase_h
#define incl_ContainerBase_h


#include <iostream>


using namespace std;


//  base class for ListBase and VectorBase

class ContainerBase
{
  public:

    ContainerBase(void) { n = 0; }

    virtual ~ContainerBase() { }

    int n;

    virtual int dim(void) { return n; }
};




#endif

