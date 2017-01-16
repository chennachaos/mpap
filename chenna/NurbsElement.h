#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#define EIGEN_SUPERLU_SUPPORT

#ifndef incl_NurbsElement_h
#define incl_NurbsElement_h


#include "List.h"
#include "Domain.h"
#include "MyString.h"
#include "DomainTypeEnum.h" // for derived elements
#include "MathVector.h"
#include "NurbsShapeFunctions.h"
#include "FunctionsMaterial.h"
#include "MpapTime.h"
//#include "FunctionsDeniz.h"
#include "NurbsSOLID.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "SolverEigen.h"

using namespace Eigen;
using namespace std;

//typedef  SparseMatrix<double>  SparseMatrixXd;
//typedef  SparseMatrix<float>   SparseMatrixXf;
//typedef  SparseMatrix<int>     SparseMatrixXi;

//typedef  DynamicSparseMatrix<double>  DynamicSparseMatrixXd;
//typedef  DynamicSparseMatrix<float>   DynamicSparseMatrixXf;
//typedef  DynamicSparseMatrix<int>     DynamicSparseMatrixXi;


class NurbsElement: public ListItem
{
  public:

    //member variables

    double  *intVar1, *intVar2, *elmDat, *matDat;

    double  uvalues[3], vvalues[3], wvalues[3], JacMultFact, elemError;

    int  nlbf, ndof, nsize, nivGP, nGP, elenum, elenum2, patchnum, counter;

    bool tracflag;

    VectorArray<double> gausspoints, gaussweights, primvar, knotsAtGPs, vals2project, resi, resi2;

    VectorArray<int>   forassy, startindex, forAssyVec1, forAssyVec2;

    ListArray<VectorArray<double> >  tracdata, stiffness_local, mass_local;

    ListArray<VectorArray<int> > forassembly, forassyKut, forassyKup, forassyKtu, forassyKtt, forassyKtp, forassyKpu, forassyKpt;

    MatrixXd  Klocal;
    VectorXd  Flocal, Flocal2;

    NurbsCURVE *curve0, *curve1;

    NurbsSURFACE *surf0, *surf1, *surf2;

    NurbsSOLID  *solid0, *solid1, *solid2;

    //member functions

    NurbsElement(void);

    virtual ~NurbsElement();

    virtual void printStiffnessMatrix();

    virtual void printForceVector();

    virtual void prepareElemData();

    virtual void printPrimVariable();

    double GetError()
    {
       return elemError;
    }

    virtual void initialiseDOFvalues()
      { cout << "   'initialiseDOFvalues' is not defined for this element!\n\n"; return; }


    virtual void initialiseKnotsAtGPs()
      { cout << "   'initialiseKnotsAtGPs' is not defined for this element!\n\n"; return; }


    virtual int calcOutput(double u1, double v1)
      { cout << "   'calcOutput' is not defined for this element!\n\n"; return 0; }

    virtual void initialiseIntVar()
      { cout << "   'initialiseIntVar' is not defined for this element!\n\n"; }

    //virtual void  AssembleElementMatrix(int, Mat, int, int)
    //{ cout << "   'AssembleElementMatrix' is not defined for this element!\n\n"; return; }

    virtual void  AssembleElementMatrix(int, SparseMatrixXd&, int, int)
    { cout << "   'AssembleElementMatrix' is not defined for this element!\n\n"; return; }

    virtual void  AssembleElementMatrix(int, SparseMatrixXd&);

    virtual void AssembleElementMatrix(int, MatrixSparseArray<double>&);
//      { cout << "   'AssembleElementMatrix' is not defined for this element!\n\n"; return; }

    virtual void AssembleElementMatrix2(int, MatrixXd&, MatrixXd& )
      { cout << "   'AssembleElementMatrix2' is not defined for this element!\n\n"; return; }

    virtual void AssembleElementMatrix3(int, double, SparseMatrixXd&);

    virtual void AssembleElementVector(bool, bool, double*, double*, int, int);

    virtual void AssembleElementVector2(bool, int, VectorXd&, double*, VectorXd&)
      { cout << "   'AssembleElementVector2' is not defined for this element!\n\n"; return; }

    virtual void diffStiffTest(double,int,int,bool)
      { cout << "   'diffStiffTest' is not defined for this element!\n\n"; return; }

    virtual int calcStiffnessMatrix(double dt)
      { cout << "   'calcStiffnessMatrix' is not defined for this element!\n\n"; return 0; }

    virtual int calcMassMatrix(int lumpInd, double dt)
      { cout << "   'calcMassMatrix' is not defined for this element!\n\n"; return 0; }

    virtual int calcInternalForces()
      { cout << "   'calcAndAssyIntForceVec' is not defined for this element!\n\n"; return 0; }

    virtual int calcLoadVector()
      { cout << "   'calcAndAssyLoadVec' is not defined for this element!\n\n"; return 0; }

    virtual int  calcStiffnessAndResidual()
      { cout << "   'calcStiffnessAndResidual' is not defined for this element!\n\n"; return 0; }

    virtual int  toComputeInfSupCondition()
      { cout << "   'toComputeInfSupCondition' is not defined for this element!\n\n"; return 0; }

    virtual void putLabel(char*, bool defFlg = false)
      { cout << "   'putLabel' is not defined for this element!\n\n"; return; }

    virtual void contourplot(int, int, double, double)
      { cout << "   'contourPlot' is not defined for this element!\n\n"; return; }

    virtual void discreteContourplot(int, int, int, int, double, double)
      { cout << "   'discreteContourplot' is not defined for this element!\n\n"; return; }

    virtual double volume(bool init = false)
      { cout << "   'volume' is not defined for this element!\n\n"; return 0.; }

    virtual void plotGaussPoints(int,bool defFlg = false)
      { cout << "  'plotGaussPoints' is not available for this element!\n\n"; return; }

    virtual bool containsPoint(double*, double*)
      { cout << "  'containsPoint' is not available for this element!\n\n"; return false; }

    virtual void projectToKnots(bool, int, int, int)
      { cout << "  'projectToNodes' is not available for this element!\n\n"; return; }

    virtual void projectIntVar(int, double*)
      { cout << "  'projectIntVar' is not available for this element!\n\n"; return; }

    virtual void projectStrain(int, int, double*)
      { cout << "  'projectStrain' is not available for this element!\n\n"; return; }

    virtual void projectStress(int, double*)
      { cout << "  'projectStress' is not available for this element!\n\n"; return; }

    virtual void createTractionDataVariable()
      { cout << "  'createTractionDataVariable' is not available for this element!\n\n"; return; }

    virtual void toPostprocess(int, int, int, SparseMatrixXd&, VectorXd& )
      { cout << "  'toPostprocess' is not available for this element!\n\n"; return; }

    virtual int applyDirichletBCs()
      { cout << "   'applyDirichletBCs' is not defined for this element!\n\n"; return 0; }

    virtual int calcError(int index)
      { cout << "   'calcError' is not defined for this element!\n\n"; return 0; }

    virtual void  computeBounds(double* val)
      { cout << "   'calcError' is not defined for this element!\n\n"; return ; }

};

#endif //incl_NurbsElement_h

