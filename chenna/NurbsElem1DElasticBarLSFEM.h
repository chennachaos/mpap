
#ifndef incl_NurbsElem1DElasticBarLSFEM_h
#define incl_NurbsElem1DElasticBarLSFEM_h


#include "NurbsElement1D.h"
#include "PropertyItem.h"


class NurbsElem1DElasticBarLSFEM: public NurbsElement
{
  public:

    double a, c, rho0, bforce;

    bool  finite;


    NurbsElem1DElasticBarLSFEM();

    //NurbsElem1DElasticBarLSFEM(NurbsShapeFns1D* shpfncs1, double gaussweight, List<PropertyItem>& );

    virtual ~NurbsElem1DElasticBarLSFEM();

    virtual void prepareElemData();

    virtual void initialiseDOFvalues();

    virtual void initialiseKnotsAtGPs();

    virtual void createTractionDataVariable();



    virtual int calcStiffnessAndResidual();

    int calcStiffnessAndResidual1();
    int calcStiffnessAndResidual2();


    virtual int calcStiffnessMatrix(double dt);

    virtual int calcMassMatrix(int lumpInd, double dt);
    
    virtual int calcInternalForces();

    virtual int calcLoadVector();

    int calcLoadVector1();
    int calcLoadVector2();

    virtual int applyDirichletBCs();

    virtual int calcOutput(double u1, double v1);

//    virtual void contourplot(int, int, double, double);


    virtual void discreteContourplot(int, int, int, int, double, double);

    virtual void projectToKnots(bool, int, int, int);

    virtual void projectStrain(int, int, double*);

    virtual void projectStress(int, double*);

    virtual void projectIntVar(int, double*);

    virtual void toPostprocess(int, int, int,  SparseMatrixXd&, VectorXd&);

    virtual void AssembleElementMatrix2(int, MatrixXd&, MatrixXd& );


};

#endif

