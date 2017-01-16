
#ifndef incl_Element2D3nodedStabIncompHighReFluid_h
#define incl_Element2D3nodedStabIncompHighReFluid_h


#include "Element2D3nodedTriangle.h"



class Element2D3nodedStabIncompHighReFluid: public Element2D3nodedTriangle
{
  public:

    Element2D3nodedStabIncompHighReFluid(void);

    virtual ~Element2D3nodedStabIncompHighReFluid();
	  
    virtual bool forDomainType(int);

    int calcStiffnessAndResidual(void);

    virtual int ndf(void) { return 3; }
    
  private:

	  
};

#endif

