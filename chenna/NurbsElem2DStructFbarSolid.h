
#ifndef incl_NurbsElem2DStructFbarSolid_h
#define incl_NurbsElem2DStructFbarSolid_h


#include "NurbsElementSolid.h"
#include "NurbsShapeFns.h"
#include "PropertyItem.h"


class NurbsElem2DStructFbarSolid: public NurbsElementSolid
{
  public:

    NurbsElem2DStructFbarSolid();

    virtual ~NurbsElem2DStructFbarSolid();
	  
    virtual int calcStiffnessAndResidual();

    int calcStiffnessAndResidual1();
    int calcStiffnessAndResidual2();
    int calcStiffnessAndResidual3();

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

    virtual void  AssembleElementMatrix(int, SparseMatrixXd&, int, int);


	  
};

#endif

