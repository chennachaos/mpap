
#ifndef incl_NurbsElem2DStokes_h
#define incl_NurbsElem2DStokes_h


#include "NurbsElementSolid.h"
#include "NurbsShapeFns.h"
#include "PropertyItem.h"


class NurbsElem2DStokes: public NurbsElementSolid
{
  public:

    NurbsElem2DStokes();

    virtual ~NurbsElem2DStokes();
	  
    virtual int calcStiffnessAndResidual();
    
    int calcStiffnessAndResidual1();
    
    int calcStiffnessAndResidual2();

    virtual int calcStiffnessMatrix(double dt);

    virtual int calcMassMatrix(int lumpInd, double dt);

    virtual int calcInternalForces();

    virtual int calcLoadVector();

    virtual void createTractionDataVariable();
    

};

#endif

