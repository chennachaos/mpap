
#ifndef incl_NurbsElement1D_h
#define incl_NurbsElement1D_h


#include "NurbsElement.h"

using namespace std;


class NurbsElement1D: public NurbsElement
{
  public:

    //member variables

  //  int startindex;

  //  NurbsCURVE *curve0, *curve1;


    //member functions
    NurbsElement1D(void);

    virtual ~NurbsElement1D();

    virtual void initialiseDOFvalues();	  
 
};

#endif //incl_NurbsElement1D_h

