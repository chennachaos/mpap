
#ifndef incl_Element1D2nodedSchroedingerCN_h
#define incl_Element1D2nodedSchroedingerCN_h


#include "Element1D2nodedLine.h"



class Element1D2nodedSchroedingerCN: public Element1D2nodedLine
{
  public:

    Element1D2nodedSchroedingerCN(void);

    virtual ~Element1D2nodedSchroedingerCN();
	  
    virtual int ndf(void) { return 2; }
    
    virtual bool forDomainType(int);

    virtual int calcStiffnessAndResidual(void);
	    
  private:

	  
};

#endif

