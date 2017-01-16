
#ifndef incl_Element1D2nodedLinearSolidALE_h
#define incl_Element1D2nodedLinearSolidALE_h


#include "Element1D2nodedLine.h"



class Element1D2nodedLinearSolidALE: public Element1D2nodedLine
{
  public:

    Element1D2nodedLinearSolidALE(void);

    virtual ~Element1D2nodedLinearSolidALE();
	  
    virtual int ndf(void) { return 3; }
    
    virtual bool forDomainType(int);

    virtual int calcStiffnessAndResidual(void);
	    
  private:

	  
};

#endif

