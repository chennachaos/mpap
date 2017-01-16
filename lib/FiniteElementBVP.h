
#ifndef incl_FiniteElementBVP_h
#define incl_FiniteElementBVP_h

#include "Mesh.h"
#include "Solver.h"


class FiniteElementBVP: public Mesh
{ 
  public:

    FiniteElementBVP(void);                      // constructor

    virtual ~FiniteElementBVP();                 // destructor

    Solver *solver;

    int localStiffnessError;
    
    virtual void readInputData(std::ifstream &, MyString &);
    
    virtual void prepareInputData(void);
    
    virtual void printInfo(void);

    virtual int  calcStiffnessAndResidual(int printRes=2, bool zeroMtx=true, bool zeroRes=true);

    virtual void setSolver(int, int *parm = NULL, bool cIO = false);

    virtual void timeUpdate(void);

    virtual void updateIterStep(void);

    virtual int  factoriseSolveAndUpdate(void);

    virtual bool converged(void);

    virtual bool diverging(double);

    virtual void addForces(void);

    virtual void elementDiffStiffTest(double,int,int,int,bool);

    virtual void elementDiffMeshDerivTest(double,int,int,int,bool);

    virtual void globalDiffStiffTest(double,int,int,bool);

    virtual void reset(void);

    virtual void printComputerTime(bool reset = true, int detailFlg = 1);

    virtual void applyCirculationCutFor2DPotentialFlow(void);

    virtual void adjustElementAssembly(int, int);

    virtual void plotSolverMatrixPattern(char*);

  private:

    double ctimFactSolvUpdt, ctimCalcStiffRes;

    bool freeCirculationFlag;

    bool checkNodalData, checkElementOutput;

};



#include "DomainInlineFunctions.h"

define_reference_cast(finiteElementBVP,FiniteElementBVP)

define_isType(isFiniteElementBVP,FINITEELEMENTBVP)
	
//define_isDerivedFrom(isDerivedFromFiniteElementBVP,FINITEELEMENTBVP)


	
#endif


