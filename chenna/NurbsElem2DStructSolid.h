
#ifndef incl_NurbsElem2DStructSolid_h
#define incl_NurbsElem2DStructSolid_h


#include "NurbsElementSolid.h"
#include "NurbsShapeFns.h"
#include "PropertyItem.h"


class NurbsElem2DStructSolid: public NurbsElementSolid
{
  public:


    NurbsElem2DStructSolid();

    virtual ~NurbsElem2DStructSolid();
	  
    virtual int calcStiffnessAndResidual();

    int calcStiffnessAndResidual1();
    int calcStiffnessAndResidual2();


    virtual int calcStiffnessMatrix(double dt);

    virtual int calcMassMatrix(int lumpInd, double dt);

    virtual int calcInternalForces();
    
    int calcInternalForces1();
    
    int calcInternalForcesAxsy();

    //virtual  void  AssembleElementMatrix(int index, Mat mtx, int start1, int start2);

    virtual void  AssembleElementMatrix(int, SparseMatrixXd&, int, int);

    virtual int calcOutput(double u1, double v1);

//    virtual void contourplot(int, int, double, double);

    virtual int calcError(int index);

    virtual void discreteContourplot(int, int, int, int, double, double);

    virtual void projectToKnots(bool, int, int, int);

    virtual void projectStrain(int, int, double*);

    virtual void projectStress(int, double*);

    virtual void projectIntVar(int, double*);

    virtual void toPostprocess(int, int, int,  SparseMatrixXd&, VectorXd&);

    virtual void  computeBounds(double* val);

};

#endif

