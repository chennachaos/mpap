
#ifndef incl_Element1D2nodedSchroedingerST_h
#define incl_Element1D2nodedSchroedingerST_h


#include "Element1D2nodedLine.h"



class Element1D2nodedSchroedingerST: public Element1D2nodedLine
{
  public:

    Element1D2nodedSchroedingerST(void);

    virtual ~Element1D2nodedSchroedingerST();
	  
    virtual int ndf(void) { return 2; }
    
    virtual bool forDomainType(int);

    virtual int calcStiffnessAndResidual(void);
	    
  private:

	  
};

#endif

