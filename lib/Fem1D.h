
#ifndef incl_Fem1D_h
#define incl_Fem1D_h

#include "FiniteElementBVP.h"



class Fem1D: public FiniteElementBVP
{ 
  public:

    Fem1D(void);                      // constructor

    virtual ~Fem1D();                 // destructor

    virtual void readInputData(std::ifstream &, MyString &);
    
    virtual void prepareInputData(void);

  private:

};



#include "DomainInlineFunctions.h"

define_reference_cast(fem1D,Fem1D)

define_isType(isFem1D,FEM1D)
	



	
#endif


