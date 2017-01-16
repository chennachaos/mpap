
#ifndef incl_NurbsElem2DStructMixed2field_h
#define incl_NurbsElem2DStructMixed2field_h


#include "NurbsElementSolid.h"
#include "NurbsShapeFns.h"
#include "PropertyItem.h"

//#define EIGEN_RUNTIME_NO_MALLOC

//Eigen::internal::set_is_malloc_allowed(true);


class NurbsElem2DStructMixed2field: public NurbsElementSolid
{
  public:

    bool calcExtraMatrices;

    MatrixXd   Kup, Kpu, Kpp;

    NurbsElem2DStructMixed2field();

    virtual ~NurbsElem2DStructMixed2field();

    virtual int calcStiffnessAndResidual();

    virtual int toComputeInfSupCondition();

    int calcStiffnessAndResidual1();
    int calcStiffnessAndResidual2();

    int calcStiffnessAndResidualAxsySS();
    int calcStiffnessAndResidualAxsyFS();

    //virtual void  AssembleElementMatrix(int, Mat, int, int);

    virtual void  AssembleElementMatrix(int, SparseMatrixXd&, int, int);

    virtual void AssembleElementMatrix(int, MatrixSparseArray<double>&);

    virtual void AssembleElementVector(bool, bool, double*, double*, int, int);

    virtual void AssembleElementMatrix2(int, MatrixXd&, MatrixXd& );

    virtual void AssembleElementVector2(bool, int, VectorXd&, double*, VectorXd&);

    virtual int calcStiffnessMatrix(double dt);

    virtual int calcMassMatrix(int lumpInd, double dt);

    virtual int calcInternalForces();

    virtual int calcOutput(double u1, double v1);

    virtual void contourplot(int, int, double, double);

    virtual int calcError(int index);

    virtual void discreteContourplot(int, int, int, int, double, double);

    virtual void projectToKnots(bool, int, int, int);

    virtual void projectStrain(int, int, double*);

    virtual void projectStress(int, double*);

    virtual void projectIntVar(int, double*);

    virtual void toPostprocess(int, int, int,  SparseMatrixXd&, VectorXd&);

};

#endif

