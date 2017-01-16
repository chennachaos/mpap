
#ifndef incl_Element2D3nodedStabCompFluid_h
#define incl_Element2D3nodedStabCompFluid_h


#include "Element2D3nodedTriangle.h"



class Element2D3nodedStabCompFluid: public Element2D3nodedTriangle
{
  public:

    Element2D3nodedStabCompFluid(void);

    virtual ~Element2D3nodedStabCompFluid();
	  
    virtual bool forDomainType(int);

    int calcStiffnessAndResidual(void);

    virtual int ndf(void) { return 3; }
    
  private:

	  
};

#endif

