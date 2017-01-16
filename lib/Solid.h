
#ifndef incl_Solid_h
#define incl_Solid_h

#include "FiniteElementBVPWNI.h"



class Solid: public FiniteElementBVPWNI
{ 
  public:
    Solid(void);                      // constructor

    virtual ~Solid();                 // destructor

    int lumpType;
    
    virtual void readInputData(std::ifstream &, MyString &);
    
    virtual void prepareInputData(void);
    
    virtual void prepareInteractions(void);
    
    virtual void printInfo(void);

    virtual void setTimeParam(void);

    virtual void timeUpdate(void);

    virtual void updateIterStep(void);

    virtual int  calcStiffnessAndResidual(int printRes=2, bool zeroMtx=true, bool zeroRes=true);

    virtual void reset(void);
    
    virtual void domainSpecificNodalDataTransfer(int);

  private:

};



#include "DomainInlineFunctions.h"

define_reference_cast(solid,Solid)

define_isType(isSolid,SOLID)
	



	
#endif


