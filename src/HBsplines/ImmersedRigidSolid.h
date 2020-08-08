#ifndef incl_ImmersedRigidSolid_h
#define incl_ImmersedRigidSolid_h


#include "headersBasic.h"
#include "ImmersedSolid.h"
#include "util.h"

using namespace std;


enum  {DOF_FIXED = -1, DOF_PRESCRIBED = -2};


class ImmersedRigidSolid : public ImmersedSolid
{
    public:

        int ndofRigidbody;
        bool  USE_SPECIFIED_PIVOT;
        MatrixXd  matM, matC, matK, Klocal, Kglobal;
        VectorXd  Flocal, rhsVec;

        vector<int>  forAssyVec, dofData;

        ImmersedRigidSolid();

        ImmersedRigidSolid(int dd);

        virtual ~ImmersedRigidSolid();

        virtual bool  isRigidBody() { return true; }
        virtual bool  isFlexibleBody() { return false; }

        virtual void  setMass(vector<double>&);

        virtual void  setDamping(vector<double>&);

        virtual void  setStiffness(vector<double>&);

        virtual  void  setNodalPositions(vector<vector<double> >&  vectemp);

        virtual void  setBoundaryConditions(vector<int>& vectemp);

        virtual void  setPrescribedMotion(vector<vector<double> >& vectemp);

        virtual void  setPivotpoint(vector<double>&);

        virtual void  setPreload(vector<double>&);

        virtual void  setInitialForcePredictor(vector<double>&);

        virtual void  setRigidBodyMotionLimits(vector<vector<double> >&);

        virtual void  printSelf();

        virtual void  initialise();

        virtual void  computeInitialAcceleration();

        virtual int  calcStiffnessAndResidual(int solver_type=1, bool zeroMtx=true, bool zeroRes=true);

        virtual int  applyBoundaryConditions();

        virtual int  applyExternalForces();

        virtual int  factoriseSolveAndUpdate();

        virtual void  solveTimeStep();

        virtual void  writeOutput();

        virtual void  postProcess(int index);

        virtual int  updatePointPositions();

        void  updatePointPositions1D();
        void  updatePointPositions2D();
        void  updatePointPositions3D();

        virtual void  updateForce();

        virtual void  updateForce(double*);

        virtual void  updateDisplacement(double*);

        virtual void  setBoundaryConditions(vector<vector<double> >& vectemp);

        virtual void setSolver(int slv = 1, int *parm = NULL, bool cIO = false);

        void  prepareMatrixPattern();

        virtual void  resetMatrixAndVector();

        virtual void calcCouplingMatrices();

        virtual int assembleGlobalMatrixAndVector(int ind1, int ind2, SparseMatrixXd& mtx, double* rhs);

        virtual int assembleGlobalMatrixAndVectorCutFEM(int ind1, int ind2, SolverPetsc* solverTemp);

        virtual void  assembleElementVector(int ind, bool flag, double* rhs);

        virtual  void  writeResult(ofstream&);

        virtual  void  readResult(ifstream&);
};






#endif
