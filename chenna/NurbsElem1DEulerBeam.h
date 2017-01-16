
#ifndef incl_NurbsElem1DEulerBeam_h
#define incl_NurbsElem1DEulerBeam_h


#include "NurbsElement1D.h"
#include "PropertyItem.h"


class  NurbsElem1DEulerBeam: public NurbsElement1D
{
  public:

    double   EI , cf,  thick,  pres;


    NurbsElem1DEulerBeam();

    virtual ~NurbsElem1DEulerBeam();

    virtual int calcStiffnessAndResidual();

    virtual int calcLoadVector();

    virtual void prepareElemData();

    virtual void initialiseDOFvalues();

    virtual void initialiseKnotsAtGPs();

    virtual int calcMassMatrix(int, double);

    virtual void AssembleElementMatrix2(int, MatrixXd&, MatrixXd& );

};

#endif

