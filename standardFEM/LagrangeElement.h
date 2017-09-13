#ifndef incl_LagrangeElement_h
#define incl_LagrangeElement_h

#include "List.h"
#include "MyString.h"
#include "MathVector.h"
#include <vector>
#include "SolverEigen.h"
#include "SolverPetsc.h"


using std::vector;
using std::cout;
using std::endl;

using Eigen::VectorXd;
using Eigen::MatrixXd;

class GeomDataLagrange;
class SolutionData;


class LagrangeElement
{
  public:

    //member variables

    bool finite, axsy, followerLoadFlag, tracflag;

    // sss=stressStrainState
    // 1.) plane stress 2.) plane strain 3.) axisymmetric

    int  elmType, matType, secType, subdomId, elmTypeNameNum;

    int  matId, finiteInt, sss, degree, npElem, ndim;

    int  nlbf, ndof, nsize, nivGP, nGP, elenum;

    double  *intVar1, *intVar2, *elmDat, *matDat;

    double  thick, elemError;

    vector<double>  resi, primvar, knotsAtGPs, vals2project, resi2;

    vector<vector<double> >  tracdata;

    vector<int>  nodeNums, forAssyVec, forAssyVec2, startindex, globalDOFnums;

    SolutionData  *SolnData;
    GeomDataLagrange  *GeomData;

    //member functions

    LagrangeElement();

    virtual ~LagrangeElement();

    int getDimension()
    { return ndim; }

    int getPolynomialDegree()
    { return degree; }

    void  setSubdomainId(int sid)
    {  subdomId = sid; return;  }

    int getSubdomainId()
    {  return  subdomId;  }

    int getNodesPerElement()
    {  return  npElem; }

    int getNdofPerNode()
    { return ndof; }

    int  getNdofPerElement()
    {  return  nsize;  }

    std::vector<int>&  getNodeNumbers()
    {  return  nodeNums; }

    std::vector<int>&  getVectorForAssembly()
    {  return  forAssyVec; }

    virtual int getElmTypeNameNum()
    {  return -1;     }

    virtual void printStiffnessMatrix();

    virtual void printForceVector();

    virtual void prepareElemData();

    virtual void prepareElemData2()
    { cout << "   'prepareElemData2()' is not defined for this element!\n\n"; return; }

    virtual void printPrimVariable();

    double getError()
    { return elemError;   }
    
    virtual void initialiseDOFvalues()
    { cout << "   'initialiseDOFvalues' is not defined for this element!\n\n"; return; }

    virtual void initialiseKnotsAtGPs()
    { cout << "   'initialiseKnotsAtGPs' is not defined for this element!\n\n"; return; }

    virtual int calcOutput(double u1, double v1)
    { cout << "   'calcOutput' is not defined for this element!\n\n"; return 0; }

    virtual void initialiseIntVar()
    { cout << "   'initialiseIntVar' is not defined for this element!\n\n"; }

    virtual void createTractionDataVariable()
    { cout << "  'createTractionDataVariable' is not available for this element!\n\n"; return; }

    virtual void diffStiffTest(double,int,int,bool)
    { cout << "   'diffStiffTest' is not defined for this element!\n\n"; return; }

    virtual int calcInternalForces()
    { cout << "   'calcAndAssyIntForceVec' is not defined for this element!\n\n"; return 0; }

    virtual int calcLoadVector()
    { cout << "   'calcAndAssyLoadVec' is not defined for this element!\n\n"; return 0; }

    virtual int  calcStiffnessAndResidual()
    { cout << "   'calcStiffnessAndResidual' is not defined for this element!\n\n"; return 0; }

    virtual int  calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal)
    { cout << "   'calcStiffnessAndResidual' is not defined for this element!\n\n"; return 0; }

    virtual int applyDirichletBCs()
    { cout << "   'applyDirichletBCs' is not defined for this element!\n\n"; return 0; }

    virtual void  assembleElementMatrixAndVector(int, SparseMatrixXd&, double*);

    virtual void assembleElementMatrix(int, Mat);

    virtual void assembleElementMatrix(int, SparseMatrixXd&);

    virtual void assembleElementVector(bool, bool, double*, double*, int start1=0, int start2=0);
    
    virtual void assembleElementVector(bool, bool, Vec, Vec, int start1=0, int start2=0);

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

    virtual  void computeEnergy(int, int, VectorXd& )
    { cout << "  'computeEnergy' is not available for this element!\n\n"; return; }

    void computeMomentum(int, int, VectorXd&);

    int calcError(int index);

    double  computeGeomOrig(int dir, VectorXd& NN);
    double  computeGeomNew(int dir, VectorXd& NN);
    double  computeGeomCur(int dir, VectorXd& NN);

    double  computeValue(int dir, VectorXd& NN);
    double  computeValuePrev(int dir, VectorXd& NN);
    double  computeValuePrev2(int dir, VectorXd& NN);
    double  computeValuePrev3(int dir, VectorXd& NN);
    double  computeValuePrev4(int dir, VectorXd& NN);
    double  computeValueExtrap(int dir, VectorXd& NN);

    double  computeValueCur(int dir, VectorXd& NN);
    double  computeValueDot(int dir, VectorXd& NN);
    double  computeValueDotCur(int dir, VectorXd& NN);

    double  computeValueDotDot(int dir, VectorXd& NN);
    double  computeValueDotDotCur(int dir, VectorXd& NN);

    double  computeValue2(int dir, VectorXd& NN);
    double  computeValue2Prev(int dir, VectorXd& NN);
    double  computeValue2Cur(int dir, VectorXd& NN);

    double  computeForce(int dir, VectorXd& NN);
    double  computeForcePrev(int dir, VectorXd& NN);
    double  computeForceCur(int dir, VectorXd& NN);
};

#endif //incl_Lagrange_Element_h

