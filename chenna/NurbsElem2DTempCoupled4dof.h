
#ifndef incl_NurbsElem2DTempCoupled4dof_h
#define incl_NurbsElem2DTempCoupled4dof_h


#include "NurbsElementSolid.h"
#include "NurbsShapeFns.h"
#include "PropertyItem.h"


class NurbsElem2DTempCoupled4dof: public NurbsElementSolid
{
  public:

    NurbsElem2DTempCoupled4dof();

    virtual ~NurbsElem2DTempCoupled4dof();
	  
    virtual int calcStiffnessAndResidual();

    virtual int calcLoadVector();

    virtual int calcStiffnessMatrix(double dt);

    virtual int calcMassMatrix(int lumpInd, double dt);

    virtual int calcInternalForces();

    virtual int calcOutput(double u1, double v1);

    virtual void createTractionDataVariable();

    virtual int calcError(int index);

};

#endif

