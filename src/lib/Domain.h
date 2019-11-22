
#ifndef incl_Domain_h
#define incl_Domain_h

#include <fstream>

#include "List.h"
#include "MyString.h"
#include "Definitions.h"
#include "MathVector.h"
#include "ElementGeom2DEnum.h"

using namespace std;



class Domain: public ListItem
{ 
  public:

    Domain(void);

    virtual ~Domain();

    void readFile(std::ifstream &);

    int ndm, ndf, tis;

    double tol, *s, *p, *xl, *ul, td[TD_DIM], ctimSinceLastCall;

    bool  solverOK, firstIter; 

    virtual void readInputData(std::ifstream &, MyString &);

    virtual void prepareInputData(void);

    virtual void prepareInteractions(void);

    // declare all member functions of all derived domains

    virtual void prepareInputData2(void)
      { cout << "  'prepareInputData2' is not available for this domain type!\n\n"; return; }

    virtual int  calcStiffnessAndResidual(int printRes=2, bool zeroMtx=true, bool zeroRes=true)
      { cout << "  'calcStiffnessAndResidual' is not available for this domain type!\n\n";
        return -1; }

    virtual void prepareElemProp(int, char**)
      { cout << "  'prepareElemProp' is not available for this domain type!\n\n"; return; }

    virtual void findMinMaxX(double*, double*, bool defm = false)
      { cout << "  'findMinMaxX' is not available for this domain type!\n\n"; return; }

    virtual void findMinMaxXBB(bool, bool, bool, bool, bool, bool, bool, double*, double*)
      { cout << "  'findMinMaxXBB' is not available for this domain type!\n\n"; return; }

    virtual void findMinMaxU(int, double&, double&)
      { cout << "  'findMinMaxU' is not available for this domain type!\n\n"; return; }

    virtual void plotMesh(bool, bool)
      { cout << "  'plotMesh' is not available for this domain type!\n\n"; return; }

    virtual void paintElemGrp(int,bool defFlg = false)
      { cout << "  'paintElemGrp' is not available for this domain type!\n\n"; return; }

    virtual void plotNodes(int, int, bool defFlg = false)
      { cout << "  'plotNodes' is not available for this domain type!\n\n"; return; }

    virtual void plotElemNum(bool defFlg = false)
      { cout << "  'plotElemNum' is not available for this domain type!\n\n"; return; }

    virtual void plotBoun(bool defFlg = false)
      { cout << "  'plotBoun' is not available for this domain type!\n\n"; return; }

    virtual void plotFixed(bool defFlg = false)
      { cout << "  'plotFixed' is not available for this domain type!\n\n"; return; }

    virtual void plotLoad(double, bool defFlg = false)
      { cout << "  'plotLoad' is not available for this domain type!\n\n"; return; }

    virtual void plotReac(double, bool defFlg = false)
      { cout << "  'plotReac' is not available for this domain type!\n\n"; return; }

    virtual void printInfo(void)
      { cout << "  'printInfo' is not available for this domain type!\n\n"; return; }

    virtual void setSolver(int, int *parm = NULL, bool cIO = false)
      { cout << "  'setSolver' is not available for this domain type!\n\n"; return; }

    virtual int  factoriseSolveAndUpdate(void)
      { cout << "  'factoriseSolveAndUpdate' is not available for this domain type!\n\n"; 
      return -1; }

    virtual void updateIterStep(void)
      { cout << "  'updateIterStep' is not available for this domain type!\n\n"; return; }

    virtual void setTimeParam(void)
      { cout << "  'setTimeParam' is not available for this domain type!\n\n"; return; }

    virtual void timeUpdate(void)
      { cout << "  'timeUpdate' is not available for this domain type!\n\n"; return; }

    virtual bool converged(void)
      { cout << "  'converged' is not available for this domain type!\n\n"; return false; }

    virtual bool diverging(double)
      { cout << "  'diverging' is not available for this domain type!\n\n"; return false; }

    virtual void addForces(void)
      { cout << "  'addForces' is not available for this domain type!\n\n"; return; }

    virtual void printNodalData(int)
      { cout << "  'printNodalData' is not available for this domain type!\n\n"; return; }

    virtual void contourPlot(int, int, int, bool, double&, double&, bool)
      { cout << "  'contourPlot' is not available for this domain type!\n\n"; return; }

    virtual void projectToNodes(MyString&, int, int)
      { cout << "  'projectToNodes' is not available for this domain type!\n\n"; return; }

    virtual void writeNodalData(void)
      { cout << "  'writeNodalData' is not available for this domain type!\n\n"; return; }

    virtual void plotGaussPoints(int,bool defFlg = false)
      { cout << "  'plotGaussPoints' is not available for this domain type!\n\n"; return; }

    virtual void plotU1D(int)
      { cout << "  'plotU1D' is not available for this domain type!\n\n"; return; }

    virtual void elementDiffStiffTest(double,int,int,int,bool)
      { cout << "  'elementDiffStiffTest' is not available for this domain type!\n\n";
        return; }

    virtual void elementDiffMeshDerivTest(double,int,int,int,bool)
      { cout << "  'elementDiffMeshDerivTest' is not available for this domain type!\n\n";
        return; }

    virtual void elementDiffStiffTestMesh(double,int,int,int,bool)
      { cout << "  'elementDiffStiffTestMesh' is not available for this domain type!\n\n";
        return; }

    virtual void globalDiffStiffTest(double,int,int,bool)
      { cout << "  'globalDiffStiffTest' is not available for this domain type!\n\n";
        return; }

    virtual void eliminateDiffTest(double,double,double,int,int,bool)
      { cout << "  'eliminateDiffTest' is not available for this domain type!\n\n";
        return; }

    virtual void updateUDepIncrements(void)
      { cout << "  'updateDependentDisplacements' is not available for this domain type!\n\n";
        return; }

    virtual void doForBending(bool, bool,bool,bool,bool,bool,bool,bool,int)
      { cout << "  'doForBending' is not available for this domain type!\n\n"; return; }

    virtual void doForSection(bool, bool, bool, bool, bool, bool, bool, bool, int, double, int)
      { cout << "  'doForSection' is not available for this domain type!\n\n"; return; }

    virtual void reset(void)
      { cout << "  'reset' is not available for this domain type!\n\n"; return; }

    virtual bool isALE(bool flag = false)
      { cout << "  'isALE' is not available for this domain type!\n\n"; return false; }

    virtual void plotInterfaceNodes(int, int, int, bool defFlg = false)
      { cout << "  'plotInterfaceNodes' is not available for this domain type!\n\n"; return; }

    virtual int  updateMesh(int, bool printRes = true)
      { cout << "  'updateMesh' is not available for this domain type!\n\n"; return -1; }

    virtual void getDisplacements(double omega = 1.)
      { cout << "  'getDisplacements' is not available for this domain type!\n\n"; return; }

    virtual void transferForces(void)
      { cout << "  'transferForces' is not available for this domain type!\n\n"; return; }

    virtual void whichSurfaces(Vector<int>&, Vector<int>&)
      { cout << "  'whichSurfaces' is not available for this domain type!\n\n"; return; }

    virtual void whichSplines(Vector<int>&, Vector<int>&)
      { cout << "  'whichSplines' is not available for this domain type!\n\n"; return; }

    virtual void remeshElemGroups(bool showFlg = false)
      { cout << "  'remeshElemGroups' is not available for this domain type!\n\n"; return; }

    virtual void calcElemSizeOpt(double, double, double, bool)
      { cout << "  'calcElemSizeOpt' is not available for this domain type!\n\n"; return; }

    virtual void calcElemSizeCurr(void)
      { cout << "  'calcElemSizeCurr' is not available for this domain type!\n\n"; return; }

    virtual void getInputElemSizeOpt(void)
      { cout << "  'getInputElemSizeOpt' is not available for this domain type!\n\n"; return; }

    virtual double getElemSizeOptInternal(double*)
      { cout << "  'getElemSizeOptInternal' is not available for this domain type!\n\n";
        return 0.; }

    virtual void smoothElemSizeOptDist(int, int, double, double)
      { cout << "  'smoothElemSizeOptDist' is not available for this domain type!\n\n"; return; }

    virtual void plotGeom(int, bool, int, bool, int*)
      { cout << "  'plotGeom' is not available for this domain type!\n\n"; return; }

    virtual void plotGeometry(unsigned int)
      { cout << "  'plotGeometry' is not available for this domain type!\n\n"; return; }

    virtual void findGeomMinMaxX(double*, double*)
      { cout << "  'findGeomMinMaxX' is not available for this domain type!\n\n"; return; }

    virtual bool elemSizeRatioOK(double)
      { cout << "  'elemSizeRatioOK' is not available for this domain type!\n\n"; return false; }

    virtual void domainSpecificNodalDataTransfer(int)
      { cout << "  'domainSpecificNodalDataTransfer' is not available for this domain type!\n\n"; 
        return; }

    virtual void transferNodalDataTest(double *)
      { cout << "  'transferNodalDataTest' is not available for this domain type!\n\n"; return; }

    virtual void writeMeshToFile(MyString &, bool defm = false)
      { cout << "  'writeMeshToFile' is not available for this domain type!\n\n"; return; }

    virtual void eliminate(bool, bool)
      { cout << "  'eliminate' is not available for this domain type!\n\n"; return; }

    virtual void prepareForExternalSolver(void*, double*, bool)
      { cout << "  'prepareForExternalSolver' is not available for this domain type!\n\n"; return; }

    virtual void setLagrangianFreeNodes(VectorArray<int> &)
      { cout << "  'setLagrangianFreeNodes' is not available for this domain type!\n\n"; return; }

    virtual void initFlightSimulatorDisplay(bool*)
      { cout << "  'initFlightSimulatorDisplay' is not available for this domain type!\n\n"; 
        return; }

    virtual void printComputerTime(bool reset = true, int detailFlg = 1);

    virtual void startComputerTime(void);

    virtual void interactiveNodeSelection(int, int, bool deselectFlag = false)
      { cout << "  'interactiveNodeSelection' is not available for this domain type!\n\n"; 
        return; }

    virtual void simpleNodeSelection(void)
      { cout << "  'simpleNodeSelection' is not available for this domain type!\n\n"; return; }

    virtual void doForVortexSheet(void)
      { cout << "  'doForVotexSheet' is not available for this domain type!\n\n"; return; }

    virtual void doForLiftingLine(void)
      { cout << "  'doForLiftingLine' is not available for this domain type!\n\n"; return; }

    virtual void doForFlexibleWing(int, bool, double *, MyString &)
      { cout << "  'doForFlexibleWing' is not available for this domain type!\n\n"; return; }

    virtual void applyCirculationCutFor2DPotentialFlow(void)
      { cout << "  'applyCirculationCutFor2DPotentialFlow' is not available for this domain type!\n\n"; return; }

    virtual void adjustElementAssembly(int, int)
      { cout << "  'adjustElementAssembly' is not available for this domain type!\n\n"; return; }

    virtual void plotSolverMatrixPattern(char*)
      { cout << "  'plotSolverMatrixPattern' is not available for this domain type!\n\n"; return; }

    virtual void solve(int, bool)
      { cout << "  'solve' is not available for this domain type!\n\n"; return; }

    virtual int doIfThisIsNeumannProb(Domain*)
      { cout << "  'doIfThisIsNeumannProb' is not available for this domain type!\n\n"; return -1; }


//------Wulf's stuff under development--------------------------------------------------------------


    virtual void strainToBoundaryDisplacement(double *)
      { cout << "  'strainToBoundaryDisplacement' is not available for this domain type!\n\n";
        return; }

    virtual void getStressFromReactions(double *)
      { cout << "  'getStressFromReactions' is not available for this domain type!\n\n"; return; }


//-----Deniz's stuff--------------------------------------------------------------------------------

    virtual void microMirror(int)
      { cout << "  'microMirror' is not available for this domain type!\n\n"; return ; }

    virtual void microGID(void)
      { cout << "  'microGID' is not available for this domain type!\n\n"; return ; }

    virtual void microDiff(double,int,int,bool,double*, int)
      { cout << "  'microDiff' is not available for this domain type!\n\n"; return ; }

    virtual void microDisp(void)
      { cout << "  'microDisp' is not available for this domain type!\n\n"; return ; }

    virtual int  microError(void)
      { cout << "  'microError' is not available for this domain type!\n\n"; return -1; }

    virtual void microRestore(void)
     { cout << "  'microRestore' is not available for this domain type!\n\n"; return ; }

    virtual int microSolve2D(double*, double*, double*, double*, bool)
      { cout << "  'microSolve2D' is not available for this domain type!\n\n"; return -1; }

    virtual int microSolve3D(double*, double*, double*, double*, bool)
      { cout << "  'microSolve3D' is not available for this domain type!\n\n"; return -1; }

    virtual void applyBoundaryDisplacement(double *, bool, bool)
      { cout << "  'strainToBoundaryDisplacement' is not available for this domain type!\n\n";
        return; }

    virtual void updatevariables(void)
      { cout << "  'updatevariables' is not available for this domain type!\n\n"; return ; }

    virtual int checkcompatibility(void)
      { cout << "  'checkcompatibility' is not available for this domain type!\n\n"; return -1; }



//-----Martina's stuff------------------------------------------------------------------------------

    virtual void solveGS(void)
      { cout << "  'solveGS' is not available for this domain type!\n\n"; return ; }



//-----Chenna's stuff-------------------------------------------------------------------------------

    virtual int deallocatePetscObjects()
      { cout << "  'deallocatePetsc' is not available for this domain type!\n\n"; return -1; }

};




#endif




