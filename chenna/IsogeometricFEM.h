
#ifndef incl_IsogeometricFEM_h
#define incl_IsogeometricFEM_h

#include "Domain.h"
#include "MathMatrix.h"
#include "List.h"
#include "MathVector.h"
#include "PropertyItem.h"
#include "PatchGroup.h"
//#include "Plot.h"
#include "SolverEigen.h"
#include "NurbsUtilitiesSURFACE.h"
#include "NurbsSOLID.h"
#include "NurbsElem2DStructSolid.h"

#include <Eigen/Dense>


using namespace Eigen;


class IsogeometricFEM: public Domain
{
  public:

    enum {INTERNAL_FORCE=0, EXTERNAL_FORCE};

    int numnp, filecount, totnumel, ntoteqs1, ntoteqs2, ntoteqs, ntotgbf1, ntotgbf2, ntotgbf, Npatch, mixedSolverFlag, SOLVER_TYPE;

    int  countcurve, countsurf, localStiffnessError, localMassError, lumpType, totIntVar, RSYS, countsolid, TSTEP;

    bool massMatrixflag, intVarFlag, globalFirstIter, kRefineFlag, defUndefFlag, comprMtxFlg, vtkFlag, LSFEMflag;

    double  rNormPrev, rNorm,  ctimFactSolvUpdt, ctimCalcStiffRes, totalError;

    double  rhoInfty;

    bool subDivStab, eqOrder;

    VectorArray<int>   assy4r, outdpatch, outddof, outdtype, assy4F, assy4F2;

    VectorArray<double>   soln, solnPrev, solnInit, solnFull, diff;

    VectorArray<double>   outdfact, reac, ForceVec, analopts;

    ListArray<VectorArray<int> > patch, patchgrpdata, tracbc2, intfdata, ibc4;

    ListArray<VectorArray<double> >  kv, dispbc, dispbc2, tracbc, forcebc, outdparam;
    ListArray<VectorArray<double> >  uu, vv, ww, outp, secVarPrev, secVarPrev2, Values;

    MatrixFullArray<double> x;

    ListArray<ListArray<EPOINT> >  Sorig, Sdef;

    List<PropertyItem>  patchElemProp, patchMatlProp;

    List<PatchGroup> patchGrp;
    
    ListArray<NurbsCURVE> CurveListOriginal, CurveListFinal, CurveResult;

    ListArray<NurbsSURFACE>  SurfaceListOriginal, SurfaceListFinal, SurfaceResult, surfSecondVar ;

    ListArray<NurbsSOLID>  SolidListOriginal, SolidListFinal, SolidResult, solidSecondVar;

    SolverEigen *solverEigen;


    NurbsElement **elem;
    
    NurbsBASE  **NurbsBaseOriginal, **NurbsBaseFinal, **NurbsBaseResult, **NurbsBaseSecondVar;

    MatrixXd  globalK, globalM, eigen_vectors, Vmat, Bmat, Qmat;

    VectorXd  eigen_values, td;


    IsogeometricFEM();

    virtual ~IsogeometricFEM();

    virtual void readInputData(std::ifstream &, MyString &);
    
    void checkInputData1();
    void checkInputData2();

    virtual void prepareInputData();

    virtual void prepareMatrixPattern();
    
    void prepareMatrixPattern1();
    void prepareMatrixPattern2();
    void prepareMatrixPattern3();
    
    void prepareMatrixPatternPetsc();

    virtual void plotGeom(int, bool, int, bool, int*);

    virtual void preparePatchElemProp();

    virtual void preparePatchMatlProp();

    virtual void processInterfaceData();

    void ProcessDispBoundaryConditions();
    
    void ProcessDispBCsCurve();
    void ProcessDispBCsSurface();
    void ProcessDispBCsSurface2();
    void ProcessDispBCsSolid();

    void ProcessBCsConstraintVariable();

    void ProcessTractionBCs();

    void ProcessTractionBCsCurves();
    void ProcessTractionBCsSurface();
    void ProcessTractionBCsSurface2();
    void ProcessTractionBCsSolid();

    void GenerateConnectivityArrays();
    
    void GenerateConnectivityArraysCurve();
    void GenerateConnectivityArraysSurface();
    void GenerateConnectivityArraysSolid();

    void GenerateConnectivityArraysConstraintVariables();

    void GenerateConnectivityArraysConstraintVariablesCurve();
    void GenerateConnectivityArraysConstraintVariablesSurface();
    void GenerateConnectivityArraysConstraintVariablesSolid();

    virtual int calcStiffnessAndResidual(int printRes=2, bool zeroMtx=true, bool zeroRes=true);

    virtual int factoriseSolveAndUpdate();

    virtual int finalsolve();

    virtual void plotGaussPoints(int,bool defFlg = false);

    virtual void PostprocessForCurves();

    virtual void PostprocessForSurfaces();

    virtual void CalcDisplacements(int);

    virtual void CalcConstraintVariables(int);

    virtual void CalcStrainStress(int);

    virtual void CalcStrainStress1(int);

    virtual void prepareInteractions();

    virtual void printInfo();

    virtual void setDifferentFlags(int, int);

    virtual void printDifferentFlags();

    virtual void printComputerTime(bool reset = true, int detailFlg = 1);

    virtual bool converged(void);

    virtual bool diverging(double);

    virtual void setTimeParam(void);

    virtual void timeUpdate(void);

    virtual void updateIterStep(void);

    virtual void findMinMaxX(double*, double*, bool defFlg = false);

    virtual void findMinMaxResult(double*, double*, bool defFlg = false);

    virtual void printData(int, int);

    virtual void contourplotStress(int, int, bool, int, int, bool, double, double, bool);

    virtual void contourplot(int, int, bool, int, int, bool, double, double, bool);

    virtual void contourplotVTK(int, int, bool, int, int, bool, double, double, int*);

    virtual void contourplotVTK2(int, int, bool, int, int, bool, double, double, int*);

    virtual void dispContourplotVTK(int, int, bool, double, double, int*);

    virtual void discreteContourplot(int, int, int, int, bool, double, double, bool);

    virtual void setSolver(int, int *parm = NULL, bool cIO = false);

    virtual void plotSolverMatrixPattern(char*);

    virtual void printResultAtPoint(int, double, double, double);

    virtual void printDispsAtParameter(int, int, double, double);

    virtual void printStrainsAtParameter(int, int, double, double);

    virtual void printStressesAtParameter(int, int, double, double);

    virtual int calcAndAssyTangentMatrix(bool flag, double dt);

    virtual int calcAndAssyLoadVector(double fact, double dt);

    virtual int calcAndAssyInternalForceVector(double dt);

    virtual void staticAnalysis(int nropt=1, int niter=1, double dt=0.0);

    virtual void quasiStaticAnalysis(int nropt=1, int niter=1, double dt=0.0);

    virtual void elementDiffStiffTest(double ddd, int elnum, int dig, int dig2, bool gfrmt);

    virtual void copyElemInternalVariables();

    virtual void writeNodalData();

//    virtual void writeNodalData1(double);

    virtual double computeReactions(int, int);

    virtual void addExternalForces();

    virtual void reset();

    virtual void readSurfaceFromFile(MyString &fileName, int geom, int patchnum);

    virtual void writeGeomToFile(MyString &fileName, int geom, int index, int patchnum);

//    virtual void writeSurfaceToFile(MyString &fileName, int geom, int patchnum);

    virtual void projectFromElemsToKnots(bool, int, int, int);

    virtual void checkContSurfs(int, int, int, int);

    virtual void printElemInvVars(int, int);

    virtual void  ModalAnalysis();

    virtual void plotFrequencyCurve(int, double);

    virtual void plotModeShape(int, int);

    virtual int BFGSsolver(int, bool, double, double);

    double  LineSearch(double, double* prsd, double* du, double STOL);

    double  gamma1(double* du, double step);

    void  postProcessFlow(int, int, int, bool, double, double, int*);
    
    void  LMSolver(int, double, double, double);
    
    void  ConvectionSolver(int, double, double, double);

    void  computeElementErrors(int);

    //void  mapper1(vtkSmartPointer<vtkUnstructuredGrid> uGrid);


  protected:

    NurbsElement* newElement(int type);

  private:

   int numPatchGrp;


};



#include "DomainInlineFunctions.h"

define_reference_cast(isogeometricFEM,IsogeometricFEM)

define_isType(isIsogeometricFEM,ISOGEOMETRICFEM)






#endif


