
#ifndef incl_NurbsElem1DAdvectionDiffusion_h
#define incl_NurbsElem1DAdvectionDiffusion_h


//#include "NurbsElement1D.h"
#include "NurbsElement.h"
#include "PropertyItem.h"


class NurbsElem1DAdvectionDiffusion: public NurbsElement
{
  public:

    double a, mu, s, tau;


    NurbsElem1DAdvectionDiffusion();

    //NurbsElem1DAdvectionDiffusion(NurbsShapeFns1D* shpfncs1, double gaussweight, List<PropertyItem>& );

    virtual ~NurbsElem1DAdvectionDiffusion();
	  
    virtual int calcStiffnessAndResidual();

    virtual int calcLoadVector();

    virtual void prepareElemData();

    virtual void initialiseDOFvalues();

    virtual void initialiseKnotsAtGPs();

//    virtual void initialiseIntVar();


};

#endif

