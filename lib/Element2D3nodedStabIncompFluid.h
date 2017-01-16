
#ifndef incl_Element2D3nodedStabIncompFluid_h
#define incl_Element2D3nodedStabIncompFluid_h


#include "Element2D3nodedTriangle.h"



class Element2D3nodedStabIncompFluid: public Element2D3nodedTriangle
{
  public:

    Element2D3nodedStabIncompFluid(void);

    virtual ~Element2D3nodedStabIncompFluid();
	  
    virtual bool forDomainType(int);

    virtual int calcStiffnessAndResidual(void);

    virtual int calcMeshDerivatives(void);

    virtual int ndf(void) { return 3; }
    
    virtual void projectError(int);

    virtual void diffMeshDerivTest(double,int,int,bool);

  private:

	  
};

#endif

