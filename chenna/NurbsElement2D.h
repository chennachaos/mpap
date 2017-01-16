
#ifndef incl_NurbsElement2D_h
#define incl_NurbsElement2D_h


#include "NurbsElement.h"

using namespace std;


class NurbsElement2D: public NurbsElement
{
  public:

    //member variables

  //  int startindex[2];

  //  NurbsSURFACE *surf0, *surf1;


    //member functions
    NurbsElement2D(void);

    virtual ~NurbsElement2D();

    virtual void initialiseDOFvalues();

    virtual void diffStiffTest(double,int,int,bool);	  
 
};

#endif //incl_NurbsElement2D_h

