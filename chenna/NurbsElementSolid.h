
#ifndef incl_NurbsElementSolid_h
#define incl_NurbsElementSolid_h


#include "NurbsElement.h"

using namespace std;
using namespace Eigen;

class NurbsElementSolid: public NurbsElement
{
  public:

    //member variables

    double thick, rho0, bforce[2];

    bool finite, axsy, followerLoadFlag;

    // sss=stressStrainState
    // 1.) plane stress 2.) plane strain 3.) axisymmetric

    int  finiteInt, matId, sss, nGP1, nGP2 ;


    //member functions
    NurbsElementSolid(void);

    virtual ~NurbsElementSolid();

    virtual void prepareElemData();

    virtual void setnivGP();

    virtual void initialiseDOFvalues();

    virtual void initialiseIntVar();

    //virtual void setGaussPointDataId();

    virtual void initialiseKnotsAtGPs();

    virtual void plotGaussPoints(int,bool defFlg = false);

    virtual int calcLoadVector();
    
    virtual void diffStiffTest(double,int,int,bool);

    virtual void createTractionDataVariable();
    
    virtual double volume(bool init = false);

//    virtual void AssembleElementMatrix(int index, MatrixXd& stiff, MatrixXd& C);
    
//    virtual void AssembleElementVector(bool, bool, VectorXd&, double*);
/*
    virtual void discreteContourplot(int, int, int, int, double, double);

    virtual void projectToKnots(bool, int, int, int);

    virtual void projectStrain(int, int, double*);

    virtual void projectStress(int, double*);

    virtual void projectIntVar(int, double*);
*/

};

#endif //incl_NurbsElementSolid_h

