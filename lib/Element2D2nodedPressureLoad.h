
#ifndef incl_Element2D2nodedPressureLoad_h
#define incl_Element2D2nodedPressureLoad_h


#include "ElementSolid.h"
#include "Element2D2nodedLine.h"



class Element2D2nodedPressureLoad: public Element2D2nodedLine
{
  public:

    Element2D2nodedPressureLoad(void);

    virtual ~Element2D2nodedPressureLoad();
	  
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

