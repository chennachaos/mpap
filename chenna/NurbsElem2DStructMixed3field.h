
#ifndef incl_NurbsElem2DStructMixed3field_h
#define incl_NurbsElem2DStructMixed3field_h


#include "NurbsElementSolid.h"
#include "NurbsShapeFns.h"
#include "PropertyItem.h"

//#define EIGEN_RUNTIME_NO_MALLOC

//Eigen::internal::set_is_malloc_allowed(true);


class NurbsElem2DStructMixed3field: public NurbsElementSolid
{
  public:

    bool calcExtraMatrices;

    MatrixXd   Kut, Kup, Ktt, Ktp;

    NurbsElem2DStructMixed3field();

    virtual ~NurbsElem2DStructMixed3field();

    virtual int calcStiffnessAndResidual();


    int calcStiffnessAndResidual1();
    int calcStiffnessAndResidual2();

    virtual void AssembleElementMatrix(int, MatrixSparseArray<double>&);

    virtual void AssembleElementVector(bool, bool, double*, double*, int, int);

    virtual void AssembleElementMatrix2(int, MatrixXd&, MatrixXd& );

    virtual void AssembleElementVector2(bool, int, VectorXd&, double*, VectorXd&);

    virtual int calcStiffnessMatrix(double dt);

    virtual int calcMassMatrix(int lumpInd, double dt);

    virtual int calcInternalForces();

    virtual int calcOutput(double u1, double v1);

    virtual void contourplot(int, int, double, double);


    virtual void discreteContourplot(int, int, int, int, double, double);

    virtual void projectToKnots(bool, int, int, int);

    virtual void projectStrain(int, int, double*);

    virtual void projectStress(int, double*);

    virtual void projectIntVar(int, double*);

    virtual void toPostprocess(int, int, int,  SparseMatrixXd&, VectorXd&);

};

#endif

