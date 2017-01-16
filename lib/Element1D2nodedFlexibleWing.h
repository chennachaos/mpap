
#ifndef incl_Element1D2nodedFlexibleWing_h
#define incl_Element1D2nodedFlexibleWing_h


#include "Element1D2nodedLine.h"



class Element1D2nodedFlexibleWing: public Element1D2nodedLine
{
  public:

    Element1D2nodedFlexibleWing(void);

    virtual ~Element1D2nodedFlexibleWing();
	  
    virtual int ndf(void) { return 3; }
    
    virtual bool forDomainType(int);

    virtual int calcStiffnessAndResidual(void);
	    
  private:

	  
};

#endif

