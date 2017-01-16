
#ifndef incl_Element1D2nodedAdvectionDiffusion_h
#define incl_Element1D2nodedAdvectionDiffusion_h


#include "Element1D2nodedLine.h"



class Element1D2nodedAdvectionDiffusion: public Element1D2nodedLine
{
  public:

    Element1D2nodedAdvectionDiffusion(void);

    virtual ~Element1D2nodedAdvectionDiffusion();
	  
    virtual int ndf(void) { return 1; }
    
    virtual bool forDomainType(int);

    virtual int calcStiffnessAndResidual(void);
	    
  private:

	  
};

#endif

