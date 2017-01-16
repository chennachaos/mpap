
#ifndef incl_RigidBody_h
#define incl_RigidBody_h


#include "Domain.h"
#include "MathVector.h"
#include "DependentDoF.h"




class RigidBody: public Domain
{ 
  public:

    RigidBody(void);
    
    virtual ~RigidBody();

    VectorArray<int>    idu, frcTmFct;
    
    VectorArray<double> prop, x, xn, x0, u, un, u3, u4, u5, u6, frc, reac;

    Vector<int>    wrndDoF, wrndType;
    Vector<double> wrndFact;

    bool    localStiffnessError;

    double  charDiameter, rNorm, rNormPrev;

    List<DependentDoF> uDep;
    
    virtual void readInputData(std::ifstream &, MyString &);

    virtual void prepareInputData(void);
    
    virtual void prepareInteractions(void);

    virtual void setTimeParam(void);

    virtual void timeUpdate(void);

    virtual void updateIterStep(void);

    virtual int  calcStiffnessAndResidual(int printRes=2, bool zeroMtx=true, bool zeroRes=true);

    virtual int  factoriseSolveAndUpdate(void);

    virtual void setSolver(int, int *parm = NULL, bool cIO = false);

    virtual bool converged(void);
   
    virtual bool diverging(double);

    virtual void writeNodalData(void);
   
  private:


};




#include "DomainInlineFunctions.h"

define_reference_cast(rigidBody,RigidBody)

define_isType(isRigidBody,RIGIDBODY)
	



#endif




