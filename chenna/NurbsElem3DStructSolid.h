
#ifndef incl_NurbsElem3DStructSolid_h
#define incl_NurbsElem3DStructSolid_h


#include "NurbsElementSolid3D.h"
#include "PropertyItem.h"


class NurbsElem3DStructSolid: public NurbsElementSolid3D
{
  public:


    NurbsElem3DStructSolid();

    virtual ~NurbsElem3DStructSolid();
	  
    virtual int calcStiffnessAndResidual();

    virtual int calcMassMatrix(int lumpInd, double dt);

    virtual int calcOutput(double u1, double v1);

//    virtual void contourplot(int, int, double, double);

    virtual int calcInternalForces();
    
    //virtual void AssembleElementMatrix(int index, Mat mtx, int start1, int start2);
    
    virtual void AssembleElementMatrix(int index, SparseMatrixXd& mtx, int start1, int start2);

    virtual void discreteContourplot(int, int, int, int, double, double);

    virtual void projectToKnots(bool, int, int, int);

    virtual void projectStrain(int, int, double*);

    virtual void projectStress(int, double*);

    virtual void projectIntVar(int, double*);

    virtual void toPostprocess(int, int, int,  SparseMatrixXd&, VectorXd&);

};

#endif

