
#ifndef incl_ElementGeomExactBeam2D_h
#define incl_ElementGeomExactBeam2D_h


#include "LagrangeElement.h"



class  ElementGeomExactBeam2D : public LagrangeElement
{
  public:

    ElementGeomExactBeam2D();

    virtual ~ElementGeomExactBeam2D();

    virtual int getElmTypeNameNum()
    {  return 4;     }

    virtual void prepareElemData();

    void prepareElemData2();

    //virtual void initialiseDOFvalues();

    virtual int  calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal);

    virtual int calcInternalForces();

    virtual int calcLoadVector();

    virtual int calcOutput(double u1, double v1);

    virtual void discreteContourplot(int, int, int, int, double, double);

    virtual void projectToKnots(bool, int, int, int);

    virtual void projectStrain(int, int, double*);

    virtual void projectStress(int, double*);

    virtual void projectIntVar(int, double*);

    virtual void toPostprocess(int, int, int,  SparseMatrixXd&, VectorXd&);
};

#endif

