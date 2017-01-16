
#ifndef incl_NurbsElem2DAdvectionDiffusion_h
#define incl_NurbsElem2DAdvectionDiffusion_h


#include "NurbsElementSolid.h"
#include "NurbsShapeFns.h"
#include "PropertyItem.h"

using namespace std;
using namespace Eigen;


class NurbsElem2DAdvectionDiffusion: public NurbsElementSolid
{
  public:

    NurbsElem2DAdvectionDiffusion();

    virtual ~NurbsElem2DAdvectionDiffusion();
	  
    virtual int calcStiffnessAndResidual();

    virtual int calcStiffnessMatrix(double dt);

    virtual int calcInternalForces();

    virtual int calcLoadVector();
    
    virtual int calcOutput(double u1, double v1);

//    virtual void contourplot(int, int, double, double);

    virtual void discreteContourplot(int, int, int, int, double, double);

    virtual void toPostprocess(int, int, int,  SparseMatrixXd&, VectorXd&);

};

#endif

