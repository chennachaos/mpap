
#ifndef incl_Fluid_h
#define incl_Fluid_h

#include "FiniteElementBVPWNI.h"


class Fluid: public FiniteElementBVPWNI
{ 
  public:
    Fluid(void);                      // constructor

    virtual ~Fluid();                 // destructor

    virtual void readInputData(std::ifstream &, MyString &);
    
    virtual void prepareInputData(void);
    
    virtual void prepareInteractions(void);
    
    virtual void printInfo(void);

    virtual void setTimeParam(void);

    virtual void timeUpdate(void);

    virtual void updateIterStep(void);

    virtual void reset(void);

  private:

    VectorArray<int> timeIntSwap;

};



#include "DomainInlineFunctions.h"

define_reference_cast(fluid,Fluid)

define_isType(isFluid,FLUID)
	

	
#endif


