
#ifndef incl_NurbsElemKirchhoffPlate_h
#define incl_NurbsElemKirchhoffPlate_h


#include "NurbsElementSolid.h"
#include "NurbsShapeFns.h"
#include "PropertyItem.h"


class NurbsElemKirchhoffPlate: public NurbsElementSolid
{
  public:


    NurbsElemKirchhoffPlate();

    virtual ~NurbsElemKirchhoffPlate();
	  
    virtual int calcStiffnessAndResidual();

    virtual int calcStiffnessMatrix(double dt);

    virtual int calcMassMatrix(int lumpInd, double dt);

    virtual int calcInternalForces();


    virtual int calcOutput(double u1, double v1);

//    virtual void contourplot(int, int, double, double);


    virtual void discreteContourplot(int, int, int, int, double, double);

    virtual void projectToKnots(bool, int, int, int);

    virtual void projectStrain(int, int, double*);

    virtual void projectStress(int, double*);

    virtual void projectIntVar(int, double*);


};

#endif

