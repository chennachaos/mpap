
#ifndef incl_LagrangeElem3DBbarFbar_h
#define incl_LagrangeElem3DBbarFbar_h


#include "LagrangeElement.h"


class  LagrangeElem3DBbarFbar : public LagrangeElement
{
  public:

    LagrangeElem3DBbarFbar();

    virtual ~LagrangeElem3DBbarFbar();

    virtual int getElmTypeNameNum()
    {  return 21; }

    virtual  void prepareElemData();

    virtual  void prepareElemData2();

    virtual int calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal);

    virtual int calcInternalForces();

    virtual int calcLoadVector();
    
    int calcStiffnessAndResidualSS(MatrixXd& Klocal, VectorXd& Flocal);

    int calcStiffnessAndResidualFS(MatrixXd& Klocal, VectorXd& Flocal);


    virtual int calcOutput(double u1, double v1);
    
    virtual void discreteContourplot(int, int, int, int, double, double);

    virtual void projectToKnots(bool, int, int, int);

    virtual void projectStrain(int, int, double*);

    virtual void projectStress(int, double*);

    virtual void projectIntVar(int, double*);

    virtual void toPostprocess(int, int, int,  SparseMatrixXd&, VectorXd&);

    virtual  void computeEnergy(int, int, VectorXd& );
};







#endif


