
#ifndef incl_FiniteElementBVPWNI_h
#define incl_FiniteElementBVPWNI_h


#include "FiniteElementBVPWI.h"



class FiniteElementBVPWNI: public FiniteElementBVPWI
{
  public:

    FiniteElementBVPWNI(void);

    virtual ~FiniteElementBVPWNI();

    bool isNeumannProblem;

    virtual void readInputData(std::ifstream &, MyString &);

    virtual void prepareInputData(void);

    virtual void prepareInteractions(void);
    
    virtual void printInfo(void);

    virtual void solve(int, bool);

    virtual void printComputerTime(bool reset = true, int detailFlg = 1);

  private:

    double ctimSolve;
 
    bool doneSolve;

    void transformIntoNeumannProblem(void);
};





#include "DomainInlineFunctions.h"

define_reference_cast(finiteElementBVPWNI,FiniteElementBVPWNI)

define_isType(isFiniteElementBVPWNI,FINITEELEMENTBVPWNI)
	


#endif




