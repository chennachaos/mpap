
#ifndef incl_NurbsElementSolid3D_h
#define incl_NurbsElementSolid3D_h


#include "NurbsElement.h"


using namespace std;


class NurbsElementSolid3D: public NurbsElement
{
  public:

    double  rho0, bforce[3];

    bool finite, followerLoadFlag;

    int  finiteInt, matId, nGP1, nGP2, nGP3, nlbf1 ;

    NurbsElementSolid3D();

    virtual ~NurbsElementSolid3D();

    virtual void prepareElemData();

    virtual void setnivGP();

    virtual void initialiseDOFvalues();

    virtual void initialiseIntVar();

    virtual void initialiseKnotsAtGPs();

    virtual void plotGaussPoints(int,bool defFlg = false);

    virtual int calcLoadVector();
    
    virtual void diffStiffTest(double,int,int,bool);

    virtual void createTractionDataVariable();
    
    virtual double volume(bool init = false);

};

#endif //incl_NurbsElementSolid3D_h

