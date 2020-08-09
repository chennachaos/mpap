#ifndef incl_HBSplineFEM_CLASS_h
#define incl_HBSplineFEM_CLASS_h


#include "headersBasic.h"
#include "HBSplineBase.h"


class HBSplineFEM: public HBSplineBase
{
    private:

        int velDOF, presDOF, fluidDOF, IBDOF, solidDOF;

        vector<vector<int> >  DDconnE, DDconnEt, DDconnG, DDconnH, DDconnHt;

    public:

        HBSplineFEM();

        virtual ~HBSplineFEM();

        virtual  void  printData(int, int);

        virtual  void  prepareInputData();

        virtual  void  printInfo();

        virtual  void  readInputData(std::ifstream &, MyString &);

        virtual  void  plotGeom(int, bool, int, bool, int*);

        void  plotGeom1D(int, bool, int, bool, int*);
        void  plotGeom2D(int, bool, int, bool, int*);
        void  plotGeom3D(int, bool, int, bool, int*);

        virtual  void  plotGaussPoints();

        void  createImmersedBoundaryPoints();

        void  ImmersedBoundaryBodyForceLSFEM();

        void  applyInterfaceTerms2D();
        void  applyInterfaceTerms3D();

        void  ImmersedBoundaryBodyForce();
        void  ImmersedBoundaryBodyForce1D();
        void  ImmersedBoundaryBodyForce2D();
        void  ImmersedBoundaryBodyForce3D();

        virtual  void  solveSolidProblem();

        virtual  void  writeNodalData();

        void  writeImmersedSolidOutput();
        void  writeFluidOutput();

        virtual  void  writeReadResult(int, string&, int);

        ///////////////////////////////////////////////////////////
        //
        // SOLUTION PHASE member functions
        //
        ///////////////////////////////////////////////////////////

        virtual  void  setTimeParam();

        virtual  void  timeUpdate();

        virtual  void  updateIterStep();

        virtual  void  reset();

        virtual  int  prepareMatrixPattern();

        void prepareMatrixPatternFluid();
        void prepareMatrixPatternSolid();
        void prepareMatrixPatternLagrangeMultipliers();

        virtual  void prepareMatrixPatternPostProcess();

        virtual  int calcStiffnessAndResidual(int printRes=2, bool zeroMtx=true, bool zeroRes=true);

        virtual  int factoriseSolveAndUpdate();

        virtual  void  addExternalForces();

        virtual  void  applyBoundaryConditions();

        virtual  void  computeElementErrors(int);

        virtual  void  setInitialConditions();

        virtual  void  computeConditionNumber();

        virtual  void  computeTotalForce(int index);

        virtual  void  printResultAtPoint(int, double, double, double);

        virtual  void  postProcessFlow(int, int, int, bool, double, double, int*);

        void  postProcess1D(int, int, int, bool, double, double, int*);
        void  postProcess2D(int, int, int, bool, double, double, int*);
        void  postProcess3D(int, int, int, bool, double, double, int*);

        void  createPostProcessGrid1D(int, int, int, bool, double, double, int*);
        void  createPostProcessGrid2D(int, int, int, bool, double, double, int*);
        void  createPostProcessGrid3D(int, int, int, bool, double, double, int*);

        virtual  void  computeTotalBodyForce(int );

        //virtual  int  solveFluidProblem();

        virtual  int  deallocatePetscObjects();
};







#include "DomainInlineFunctions.h"

define_reference_cast(hbsplineFEM, HBSplineFEM)

define_isType(isHBSplineFEM, HBSPLINEFEM)


#endif









