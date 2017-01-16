
#ifndef incl_Element1D2nodedPipeFlowST_h
#define incl_Element1D2nodedPipeFlowST_h


#include "Element1D2nodedLine.h"



class Element1D2nodedPipeFlowST: public Element1D2nodedLine
{
  public:

    Element1D2nodedPipeFlowST(void);

    virtual ~Element1D2nodedPipeFlowST();
	  
    virtual int ndf(void) { return 3; }
    
    virtual bool forDomainType(int);

    virtual int calcStiffnessAndResidual(void);
	    
  private:

	  
};

#endif

