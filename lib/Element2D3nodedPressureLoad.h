
#ifndef incl_Element2D3nodedPressureLoad_h
#define incl_Element2D3nodedPressureLoad_h


#include "ElementSolid.h"
#include "Element2D3nodedLine.h"



class Element2D3nodedPressureLoad: public Element2D3nodedLine
{
  public:

    Element2D3nodedPressureLoad(void);

    virtual ~Element2D3nodedPressureLoad();
	  
    virtual bool forDomainType(int);

    int calcStiffnessAndResidual(void);

    virtual int nGaussPoints(void) { return 0; }
    
    virtual int ndf(void) { return 2; }

    virtual double volume(bool init = false) { return 0.; }

    virtual void initialiseIntVar(void) { return; }

    virtual int finiteStrain(void);

  private:

	  
};

#endif

