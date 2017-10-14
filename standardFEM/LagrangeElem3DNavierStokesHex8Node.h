
#ifndef incl_LagrangeElem3DNavierStokesHex8Node_h
#define incl_LagrangeElem3DNavierStokesHex8Node_h


#include "LagrangeElement.h"


class  LagrangeElem3DNavierStokesHex8Node : public LagrangeElement
{
  public:

    LagrangeElem3DNavierStokesHex8Node();

    virtual ~LagrangeElem3DNavierStokesHex8Node();

    virtual int getElmTypeNameNum()
    {  return 25; }

    virtual void prepareElemData();

    void prepareElemData2();

    virtual int calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal, bool firstIter=false);

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

