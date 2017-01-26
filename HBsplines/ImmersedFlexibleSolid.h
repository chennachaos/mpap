#ifndef incl_ImmersedFlexibleSolid_h
#define incl_ImmersedFlexibleSolid_h

#include "Solver.h"

#include "ImmersedSolid.h"

//#include "util.h"

#include "TreeNode.h"

#include "LagrangeElement.h"



class ImmersedFlexibleSolid : public ImmersedSolid
{
    public:

        vector<vector<double> >  DirichletBCs;

        vector<vector<int> >  NodeType, IEN, ID, LM;

        vector<int>  node_map_new_to_old, node_map_old_to_new;
        vector<int>  dof_map_new_to_old, dof_map_old_to_new;


        LagrangeElement  **elems;

        ImmersedFlexibleSolid();
        
        ImmersedFlexibleSolid(int dd);

        virtual ~ImmersedFlexibleSolid();

        virtual bool IsRigidBody() { return false; }
        virtual bool IsFlexibleBody() { return true; }

        void  prepareElemProp();
        void  prepareMatlProp();

        virtual void  printInfo();

        virtual void  initialise();

        virtual void  computeInitialAcceleration();

        virtual  void  SetNodalPositions(vector<vector<double> >&  vectemp);

        virtual int   calcStiffnessAndResidual(int solver_type, bool zeroMtx, bool zeroRes);

        virtual void  calcForceVector();

        virtual void  applyBoundaryConditions(int start, SparseMatrixXd& globalK, double* rhs);

        virtual void  applyExternalForces();

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

        int  factoriseSolveAndUpdate();

        virtual void calcCouplingMatrices();

        virtual int AssembleGlobalMatrixAndVector(int ind1, int ind2, SparseMatrixXd& mtx, double* rhs);

        virtual int AssembleGlobalMatrixAndVectorCutFEM(int ind1, int ind2, SparseMatrixXd& mtx, double* rhs);

        virtual void setSolver(int slv = 1, int *parm = NULL, bool cIO = false);

        virtual void  SetSolidElements(vector<vector<int> >& datatemp);

        virtual void  SetBoundaryConditions(vector<vector<double> >& vectemp);

        void  prepareMatrixPattern();
};






#endif