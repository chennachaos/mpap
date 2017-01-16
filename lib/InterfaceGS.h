
#ifndef incl_InterfaceGS_h
#define incl_InterfaceGS_h


#include "InterfaceMatch.h"


class InterfaceGS: public InterfaceMatch
{ 
  public:

    InterfaceGS(void);

    virtual ~InterfaceGS();

    virtual void readInputData(std::ifstream &, MyString &);

    virtual void prepareInputData(void);

    virtual void prepareInteractions(void);

    virtual void printInfo(void);

    virtual int  doIfThisIsNeumannProb(Domain*);

    virtual void setSolver(int, int *parm = NULL, bool cIO = false);

    virtual int  calcStiffnessAndResidual(int printRes=2, bool zeroMtx=true, bool zeroRes=true);

    virtual int  factoriseSolveAndUpdate(void);

    virtual bool converged(void);

    virtual void timeUpdate(void);

  private:

    int neumannDom, relaxMode, nMyAitken;

    double relaxParam[5];

    VectorArray<double> trac, tracPrev, mtx, dmx, rhs, alp;

    VectorArray<int> pos;

    List< VectorArray<double> > ua, Du;

    VectorArray<double> uInc, DDui, DDuj;

    int getAlphaForMyAitken(void);
};


enum { NONE, FIXED, AITKEN, MYAITKEN };



#include "DomainInlineFunctions.h"

define_reference_cast(interfaceGS,InterfaceGS)

define_isType(isInterfaceGS,INTERFACEGS)

	
#endif


