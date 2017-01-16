
#ifndef incl_NurbsElem2DStructSolidLSFEM3dof_h
#define incl_NurbsElem2DStructSolidLSFEM3dof_h


#include "NurbsElementSolid.h"
#include "NurbsShapeFns.h"
#include "PropertyItem.h"


class NurbsElem2DStructSolidLSFEM3dof: public NurbsElementSolid
{
  public:
    
    MatrixXd   Kup, Kpp;

    NurbsElem2DStructSolidLSFEM3dof();

    virtual ~NurbsElem2DStructSolidLSFEM3dof();

    virtual int calcStiffnessAndResidual();
    
    virtual int toComputeInfSupCondition();

    int calcStiffnessAndResidual1();
    int calcStiffnessAndResidual2();


    virtual int calcStiffnessMatrix(double dt);

    virtual int calcMassMatrix(int lumpInd, double dt);
    
    virtual int calcInternalForces();

    /*    
    int calcInternalForces1();
    
    int calcInternalForcesAxsy();
    */

    virtual int calcLoadVector();

    int calcLoadVector1();
    int calcLoadVector2();

    virtual int applyDirichletBCs();

    virtual int calcOutput(double u1, double v1);

//    virtual void contourplot(int, int, double, double);


    virtual void discreteContourplot(int, int, int, int, double, double);

    virtual void projectToKnots(bool, int, int, int);

    virtual void projectStrain(int, int, double*);

    virtual void projectStress(int, double*);

    virtual void projectIntVar(int, double*);

    virtual void toPostprocess(int, int, int,  SparseMatrixXd&, VectorXd&);

    virtual void createTractionDataVariable();

};

#endif

