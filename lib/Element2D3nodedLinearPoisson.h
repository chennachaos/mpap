
#ifndef incl_Element2D3nodedLinearPoisson_h
#define incl_Element2D3nodedLinearPoisson_h


#include "Element2D3nodedTriangle.h"



class Element2D3nodedLinearPoisson: public Element2D3nodedTriangle
{
  public:

    Element2D3nodedLinearPoisson(void);

    virtual ~Element2D3nodedLinearPoisson();
	  
    virtual bool forDomainType(int);

    int calcStiffnessAndResidual(void);

    virtual int ndf(void) { return 1; }

  private:

};

#endif

