
#ifndef incl_SolidALE_h
#define incl_SolidALE_h

#include "FiniteElementBVPWI.h"



class SolidALE: public FiniteElementBVPWI
{ 
  public:

    SolidALE(void);                      // constructor

    virtual ~SolidALE();                 // destructor

    virtual void readInputData(std::ifstream &, MyString &);
    
    virtual void prepareInputData(void);
    
    virtual void prepareInteractions(void);
    
    virtual void printInfo(void);

    virtual void setTimeParam(void);

    virtual void timeUpdate(void);

    virtual void updateIterStep(void);

    virtual int  calcStiffnessAndResidual(int printRes=2, bool zeroMtx=true, bool zeroRes=true);

  private:

    //Transport *transport;
};



#include "DomainInlineFunctions.h"

define_reference_cast(solidALE,SolidALE)

define_isType(isSolidALE,SOLIDALE)
	



	
#endif


