#ifndef incl_ImmersedFlexibleSolid_h
#define incl_ImmersedFlexibleSolid_h

#include "Solver.h"
#include "ImmersedSolid.h"
#include "TreeNode.h"
#include "LagrangeElement.h"



class ImmersedFlexibleSolid : public ImmersedSolid
{
    public:

        vector<vector<double> >  DirichletBCs;

        vector<vector<int> >  NodeType, IEN, ID, LM;

        vector<int>  node_map_new_to_old, node_map_old_to_new;


        LagrangeElement  **elems;

        ImmersedFlexibleSolid();

        ImmersedFlexibleSolid(int dd);

        virtual ~ImmersedFlexibleSolid();

        virtual bool isRigidBody() { return false; }
        virtual bool isFlexibleBody() { return true; }

        void  prepareElemProp();
        void  prepareMatlProp();
        
        int   setBoundaryConditions();

        virtual void  printInfo();

        virtual void  initialise();

        virtual void  computeInitialAcceleration();

        virtual void  setNodalPositions(vector<vector<double> >&  vectemp);

        virtual int   calcStiffnessAndResidual(int solver_type, bool zeroMtx, bool zeroRes);

        virtual void  calcForceVector();

        virtual int  applyBoundaryConditions(int start1, int start2, SparseMatrixXd& globalK, double* rhs);

        virtual int  applyBoundaryConditions(int start1, int start2, Mat mtxTemp, Vec rhsTemp);

        virtual int  applyExternalForces();

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

        int  factoriseSolveAndUpdate();

        virtual void calcCouplingMatrices();

        virtual int assembleGlobalMatrixAndVector(int ind1, int ind2, SparseMatrixXd& mtx, double* rhs);

        virtual int assembleGlobalMatrixAndVectorCutFEM(int ind1, int ind2, SolverPetsc* solverTemp);

        virtual void setSolver(int slv = 1, int *parm = NULL, bool cIO = false);

        virtual void  setSolidElements(vector<vector<int> >& datatemp);

        virtual void  setBoundaryConditions(vector<vector<double> >& vectemp);

        void  prepareMatrixPattern();

        void  addControlTerms();

        virtual  void  writeResult(ofstream&);

        virtual  void  readResult(ifstream&);
};






#endif
