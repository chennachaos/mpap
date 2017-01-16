
#ifndef incl_NurbsElem2DHeatTransfer_h
#define incl_NurbsElem2DHeatTransfer_h


#include "NurbsElementSolid.h"
#include "NurbsShapeFns.h"
#include "PropertyItem.h"


class NurbsElem2DHeatTransfer: public NurbsElementSolid
{
  public:

    NurbsElem2DHeatTransfer();

    virtual ~NurbsElem2DHeatTransfer();
	  
    virtual int calcStiffnessAndResidual();

    virtual int calcStiffnessMatrix(double dt);

    virtual int calcInternalForces();

    virtual int calcLoadVector();
    
    virtual int calcOutput(double u1, double v1);

//    virtual void contourplot(int, int, double, double);

    virtual void  AssembleElementMatrix(int, SparseMatrixXd&, int, int);

    virtual void discreteContourplot(int, int, int, int, double, double);

    virtual void toPostprocess(int, int, int,  SparseMatrixXd&, VectorXd&);

    virtual void projectToKnots(bool, int, int, int);

    virtual void projectStrain(int, int, double*);

    virtual int calcError(int index);

};

#endif

