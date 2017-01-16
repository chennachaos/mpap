#ifndef incl_LagrangeElement_h
#define incl_LagrangeElement_h

//#include "headersBasic.h"

#include "List.h"
//#include "Domain.h"
#include "MyString.h"
//#include "DomainTypeEnum.h" // for derived elements
#include "MathVector.h"

#include <vector>

//#include "FunctionsMaterial.h"
//#include "MpapTime.h"
//#include "FunctionsDeniz.h"

#include "SolverEigen.h"
#include "SolverPetsc.h"

#include <Eigen/Dense>

//using namespace Eigen;
//using namespace std;

using std::vector;
using std::cout;
using std::endl;

using Eigen::VectorXd;
using Eigen::MatrixXd;

class GeomDataLagrange;
class SolutionData;
//class SolverEigen;
//class SolverPetsc;


class LagrangeElement
{
  public:

    //member variables

    bool finite, axsy, followerLoadFlag, tracflag;

    // sss=stressStrainState
    // 1.) plane stress 2.) plane strain 3.) axisymmetric

    int  elmType, matType, secType, subdomId;

    int  matId, finiteInt, sss, degree, npElem, ndim;

    int  nlbf, ndof, nsize, nivGP, nGP, elenum;

    double  *intVar1, *intVar2, *elmDat, *matDat;

    double  thick, elemError;

    vector<double>  resi, primvar, knotsAtGPs, vals2project, resi2;

    vector<vector<double> >  tracdata;

    vector<int>  nodeNums, forAssyVec, forAssyVec2, startindex, globalDOFnums;

    //VectorXd  Flocal, Flocal2;
    //MatrixXd  Klocal;

    SolutionData  *SolnData;
    GeomDataLagrange  *GeomData;

    //member functions

    LagrangeElement(void);

    virtual ~LagrangeElement();

    int GetDOFPerNode()
    { return ndof; }

    int GetPolynomicalDegree()
    { return degree; }

    int GetDimension()
    { return ndim; }

    void  set_subdomain_id(int sid)
    {  subdomId = sid; return;  }

    int get_subdomain_id()
    {  return  subdomId;  }

    int get_nodes_per_element()
    {  return  npElem; }

    //int get_ndof_per_node()
    //{  return  ndof;  }

    int  get_ndof_per_element()
    {  return  nsize;  }

    std::vector<int>&  get_node_numbers()
    {  return  nodeNums; }

    std::vector<int>&  get_vector_for_assembly()
    {  return  forAssyVec; }


    virtual void printStiffnessMatrix();

    virtual void printForceVector();

    virtual void prepareElemData();

    virtual void prepareElemData2()
      { cout << "   'prepareElemData2()' is not defined for this element!\n\n"; return; }

    virtual void printPrimVariable();

    double GetError()
     { return elemError;   }
    
    virtual void  resetMatrixAndVector()
      { cout << "   'resetMatrixAndVector' is not defined for this element!\n\n"; return; }

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


    virtual void  AssembleElementMatrixAndVector(int, SparseMatrixXd&, double*);
      //{ cout << "   'AssembleMatrixAndVector' is not defined for this element!\n\n"; return; }

    virtual void AssembleElementMatrix(int, Mat);
      //{ cout << "   'AssembleElementMatrix' is not defined for this element!\n\n"; return; }

    virtual void AssembleElementMatrix(int, SparseMatrixXd&);

    virtual void AssembleElementVector(bool, bool, double*, double*, int start1=0, int start2=0);
    //{ cout << "   'AssembleElementVector' is not defined for this element!\n\n"; return; }
    
    virtual void AssembleElementVector(bool, bool, Vec, Vec, int start1=0, int start2=0);

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
    //{ cout << "  'computeMomentum' is not available for this element!\n\n"; return; }

    //virtual void toPostprocess(int, int, int, SparseMatrixXd&, VectorXd& )
      //{ cout << "  'toPostprocess' is not available for this element!\n\n"; return; }

    int calcError(int index);
    //{ cout << "   'calcError' is not defined for this element!\n\n"; return 0; }


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

