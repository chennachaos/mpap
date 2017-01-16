
#ifndef incl_NurbsElem1DElasticBar_h
#define incl_NurbsElem1DElasticBar_h


#include "NurbsElement1D.h"
#include "PropertyItem.h"


class NurbsElem1DElasticBar: public NurbsElement1D
{
  public:

    double   a, c, c0, rho0, bforce;

    bool  finite;


    NurbsElem1DElasticBar();

    virtual ~NurbsElem1DElasticBar();

    virtual int calcStiffnessAndResidual();

    virtual int calcLoadVector();

    virtual void prepareElemData();

    virtual void initialiseDOFvalues();

    virtual void initialiseKnotsAtGPs();

    virtual int calcMassMatrix(int, double);

    virtual void AssembleElementMatrix2(int, MatrixXd&, MatrixXd& );
};

#endif

