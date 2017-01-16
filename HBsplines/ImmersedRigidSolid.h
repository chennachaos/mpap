#ifndef incl_ImmersedRigidSolid_h
#define incl_ImmersedRigidSolid_h


#include "headersBasic.h"
#include "ImmersedSolid.h"
#include "util.h"

using namespace std;


class ImmersedRigidSolid : public ImmersedSolid
{
    public:

        int size;
        MatrixXd  matM, matC, matK, Klocal, Kglobal;
        VectorXd  Flocal, rhsVec;

        vector<int>  forAssyVec, dofData;

        ImmersedRigidSolid();

        ImmersedRigidSolid(int dd);

        virtual ~ImmersedRigidSolid();

        virtual bool  IsRigidBody() { return true; }
        virtual bool  IsFlexibleBody() { return false; }

        virtual void  SetMass(vector<double>&);

        virtual void  SetDamping(vector<double>&);

        virtual void  SetStiffness(vector<double>&);

        virtual  void  SetNodalPositions(vector<vector<double> >&  vectemp);

        virtual void  SetBoundaryConditions(vector<int>& vectemp);

        virtual void  SetPrescibedMotion(vector<int>& vectemp);

        virtual void  printSelf();

        virtual void  initialise();

        virtual void  computeInitialAcceleration();

        virtual int  calcStiffnessAndResidual(int solver_type=1, bool zeroMtx=true, bool zeroRes=true);

        virtual void  applyBoundaryConditions();

        virtual void  applyExternalForces();

        virtual int  factoriseSolveAndUpdate();

        virtual void  SolveTimeStep();

        virtual void  writeOutput();

        virtual void  postProcess(int index);

        virtual int  updatePointPositions();

        void  updatePointPositions1D();
        void  updatePointPositions2D();
        void  updatePointPositions3D();

        virtual void  updateForce();

        virtual void  updateForce(double*);

        virtual void  updateDisplacement(double*);

        virtual void  SetBoundaryConditions(vector<vector<double> >& vectemp);

        virtual void setSolver(int slv = 1, int *parm = NULL, bool cIO = false);

        void  prepareMatrixPattern();

        virtual void  resetMatrixAndVector();

        virtual void calcCouplingMatrices();

        virtual int AssembleGlobalMatrixAndVector(int ind1, int ind2, SparseMatrixXd& mtx, double* rhs);

        virtual int AssembleGlobalMatrixAndVectorCutFEM(int ind1, int ind2, SparseMatrixXd& mtx, double* rhs);

        virtual void  AssembleElementVector(int ind, bool flag, double* rhs);
};






#endif
