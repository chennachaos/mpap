#ifndef incl_ImmersedSolid_h
#define incl_ImmersedSolid_h

#include "SolutionData.h"
#include "GeomDataLagrange.h"
#include "AABB.h"
#include "Ray.h"
#include "myPoly.h"
#include "headersVTK.h"
#include  "SolverEigen.h"

#include <vector>
#include <iostream>

using std::vector;
using std::cout;
using std::endl;

using namespace myGeom;

class  ImmersedIntegrationElement;


enum  {BC_ENFORCE_TYPE_LAGRANGE=0, BC_ENFORCE_TYPE_PENALTY};



class ImmersedSolid
{
    public:
    
        static  int  count;

        int  id, DIM, type, localStiffnessError, filecount, iterCount, tis, ndof, solverOK ;
        int  nElem, nImmInt, totalDOF, npElem, nNode, nsize;
        int  matId, elemId, BC_ENFORCE_TYPE;

        double  rNorm, rNormPrev, tol, dt, rho, beta, ctimCalcStiffRes, ctimFactSolvUpdt, PENALTY, NitscheFact;

        bool firstIter, STAGGERED, isNitsche, PRESC_MOTION;

        vector<int>  assy4r, PrescMotionTimeFuncs;

        vector<vector<int> >  forAssyMat, forAssyCoupledHorz, forAssyCoupledVert;
        vector<vector<int> >  OutputData;

        VectorXd  soln, totalForce;
        myPoint  centroid;

        MatrixXd  Khorz, Kvert;

        SolutionData  SolnData;
        GeomDataLagrange  GeomData;
        
        AABB  bbox;

        Ray ray1;

        vector<ImmersedIntegrationElement*>  ImmIntgElems;
        
        vector<myPoly*>  ImmersedFaces;

        SolverEigen *solver;

        vtkSmartPointer<vtkUnstructuredGrid>  uGrid;

        vtkSmartPointer<vtkPolyData>  polyDataVTK;

        vtkSmartPointer<vtkSelectEnclosedPoints> selectEnclosedPoints;
        vtkSmartPointer<vtkSelectEnclosedPoints> selectEnclosedPoints2;

        ImmersedSolid();

        virtual ~ImmersedSolid();

        int GetID()
        {  return id; }

        void  SetDimension(int dd)
        { DIM = dd; }

        void  SetNNodes(int nn)
        {  nNode = nn; }

        int  GetNumNodes()
        { return nNode; } 

        void  SetNelem(int nn)
        {  nElem = nn; }

        int  GetNelem()
        {  return nElem; }

        void  SetNodePerElem(int nn)
        {  npElem = nn; }

        int GetTotalDOF()
        { return  totalDOF; }
        
        void  SetTolerance(double tt)
        { tol = tt; }

        void  SetPenaltyParameter(double tt)
        { PENALTY = tt; }

        double  GetPenaltyParameter()
        { return PENALTY; }

        void  SetNitscheFact(double tt)
        { NitscheFact = tt; }

        double  GetNitscheFact()
        { return NitscheFact; }

        void  SetNitscheFlag(bool tt)
        { isNitsche = tt; }

        bool  GetNitscheFlag()
        { return isNitsche; }

        void SetTimeIncrementType(int ttt)
        {  SolnData.SetTimeIncrementType(ttt); }

        void SetRho(double ttt)
        {  SolnData.SetRho(ttt); }

        void SetInitDisplacement(int ind, int dir, double ttt)
        {  SolnData.var1[ind*ndof+dir] = ttt; }
        
        void SetInitVelocity(int ind, int dir, double ttt)
        {  SolnData.var1Dot[ind*ndof+dir] = ttt; }

        void SetInitAcceleration(int ind, int dir, double ttt)
        {  SolnData.var1DotDot[ind*ndof+dir] = ttt; }

        double GetDisplacement(int ind, int dir)
        {  return  SolnData.var1[ind*ndof+dir]; }

        double GetDisplacementCur(int ind, int dir)
        {  return  SolnData.var1Cur[ind*ndof+dir]; }

        double GetVelocity(int ind, int dir)
        {  return  SolnData.var1Dot[ind*ndof+dir]; }

        double GetVelocityCur(int ind, int dir)
        {  return  SolnData.var1DotCur[ind*ndof+dir]; }

        double GetAcceleration(int ind, int dir)
        {  return  SolnData.var1DotDot[ind*ndof+dir]; }

        double GetAccelerationCur(int ind, int dir)
        {  return  SolnData.var1DotDotCur[ind*ndof+dir]; }

        double GetForce(int ind, int dir)
        {  return  SolnData.force[ind*ndof+dir]; }

        double GetForceCur(int ind, int dir)
        {  return  SolnData.forceCur[ind*ndof+dir]; }

        void  SetBoundaryConditionType(int tt)
        {
          if(tt)
            BC_ENFORCE_TYPE = BC_ENFORCE_TYPE_LAGRANGE;
          else
            BC_ENFORCE_TYPE = BC_ENFORCE_TYPE_PENALTY;
        }

        bool IsBoundaryConditionTypeLagrange()
        {
          return (BC_ENFORCE_TYPE == BC_ENFORCE_TYPE_LAGRANGE);
        }

        void  SetMaterialID(int tt)
        { matId = tt; }

        int  GetMaterialID()
        { return matId; }

        void  SetElementID(int tt)
        { elemId = tt; }

        int  GetElementID()
        { return elemId; }

        virtual bool IsRigidBody()
        { cout << "   'IsRigidBody()' is not defined for this Solid!\n\n"; return true; }

        virtual bool IsFlexibleBody()
        { cout << "   'IsFlexibleBody()' is not defined for this Solid!\n\n"; return true; }

        virtual  void  SetNodalPositions(vector<vector<double> >&  vectemp) = 0;
        //{ cout << "   'SetNodalPositions' is not defined for this Solid!\n\n"; return; }

        virtual void  SetSolidElements(vector<vector<int> >& datatemp)
        { cout << "   'SetSolidElements' is not defined for this Solid!\n\n"; return; }

        virtual void  SetImmersedElemActiveFlag(vector<int>& datatemp);
        //{ cout << "   'SetImmersedElemActiveFlag' is not defined for this Solid!\n\n"; return; }

        virtual void  SetImmersedIntegrationElements(vector<vector<int> >& datatemp);

        virtual  void  SetImmersedFaces();
        
        virtual  void  UpdateImmersedFaces();

        void  adjustBoundaryPoints(double* minVal, double* maxVal);

        virtual void  SetDataForOutput(vector<vector<int> >& vectemp);

        virtual void  SetBoundaryConditions(vector<vector<double> >& vectemp)
        { cout << "   'SetBoundaryConditions' is not defined for this Solid!\n\n"; return; }

        virtual void  printSelf()
        { cout << "   'printSelf()' is not defined for this Solid!\n\n"; return; }

        virtual void  initialise()
        { cout << "   'initialise()' is not defined for this Solid!\n\n"; return; }

        virtual void  SetMass(vector<double>&)
        { cout << "   'SetMass()' is not defined for this Solid!\n\n"; return; }

        virtual void  SetDamping(vector<double>&)
        { cout << "   'SetDamping()' is not defined for this Solid!\n\n"; return; }

        virtual void  SetStiffness(vector<double>&)
        { cout << "   'SetStiffness()' is not defined for this Solid!\n\n"; return; }

        virtual void  SetBoundaryConditions(vector<int>& vectemp)
        { cout << "   'SetBoundaryConditions' is not defined for this Solid!\n\n"; return; }

        virtual void  SetPrescibedMotion(vector<int>& vectemp)
        { cout << "   'SetPrescibedMotion' is not defined for this Solid!\n\n"; return; }

        virtual void SetSolver(int slv = 1, int *parm = NULL, bool cIO = false)
        { cout << "   'SetSolver()' is not defined for this Solid!\n\n"; return; }

        virtual void  computeInitialAcceleration()
        { cout << "   'computeInitialAcceleration()' is not defined for this Solid!\n\n"; return; }

        virtual int  calcStiffnessAndResidual(int solver_type=1, bool zeroMtx=true, bool zeroRes=true)
        { cout << "   'calcStiffnessAndResidual()' is not defined for this Solid!\n\n"; return 0; }
        
        virtual void  applyBoundaryConditions(int start1, SparseMatrixXd& globalK, double* rhs)
        { cout << "   'applyBoundaryConditions()' is not defined for this Solid!\n\n"; return; }
        
        virtual void  applyExternalForces()
        { cout << "   'applyExternalForces()' is not defined for this Solid!\n\n"; return; }

        virtual int  factoriseSolveAndUpdate()
        { cout << "   'factoriseSolveAndUpdate()' is not defined for this Solid!\n\n"; return 0; }

        virtual void  setTimeParam();

        virtual bool  diverging(double factor);

        virtual bool  converged();

        virtual void  SolveTimeStep()
        { cout << "   'SolveTimeStep()' is not defined for this Solid!\n\n"; return; }

        virtual void  timeUpdate();

        virtual void  updateIterStep();

        virtual void  reset();

        void  computeCentroid(int index);

        myPoint&  GetCentroid(int index)
        {
          computeCentroid(index);

          return  centroid;
        }

        void  computeAABB(int index);

        virtual  int  within(myPoint& pt);
        
        virtual  int  doAABBintersect(AABB&  bb2);

        virtual  int  doIntersect2D(AABB&  bb2, bool  flag, vector<int>&  vecTemp, vector<myPoint>&  ptOut);

        virtual  int  doIntersect2Dfor3D(int sideTemp, double coord3, AABB& bbTemp, bool flag, vector<int>&  vecTemp,  vector<myPoint>& ptOut);

        virtual  int  doIntersect3D(AABB&  bb2, bool  flag, vector<int>&  vecTemp, vector<myPoint>&  ptOut);

        double  distanceFromPoint(myPoint&  pt);

        double  distanceFromPoint(double xx=0.0, double yy=0.0, double zz=0.0);

        void  computeTotalForce();

        virtual void calcCouplingMatrices()
        { cout << "   'calcCouplingMatrices()' is not defined for this Solid!\n\n"; return; }

        virtual int AssembleGlobalMatrixAndVector(int ind1, int ind2, SparseMatrixXd& mtx, double* rhs)
        { cout << "   'AssembleGlobalMatrixAndVector()' is not defined for this Solid!\n\n"; return 0; }

        virtual int AssembleGlobalMatrixAndVectorCutFEM(int ind1, int ind2, SparseMatrixXd& mtx, double* rhs)
        { cout << "   'AssembleGlobalMatrixAndVectorCutFEM()' is not defined for this Solid!\n\n"; return 0; }

        virtual void  writeOutput()
        { cout << "   'writeOutput()' is not defined for this Solid!\n\n"; return; }

        virtual void  postProcess(int index)
        { cout << "   'postProcess()' is not defined for this Solid!\n\n"; return; }

        virtual int  updatePointPositions()
        { cout << "   'updatePointPositions()' is not defined for this Solid!\n\n"; return 0; }

        virtual void  updateForce() // routine for rigid-bodies with FDM
        { cout << "   'updateForce()' is not defined for this Solid!\n\n"; return; }

        virtual void  updateForce(double*) // routine for rigid-bodies with CutFEM
        { cout << "   'updateForce(double*)' is not defined for this Solid!\n\n"; return; }

        virtual void  updateDisplacement(double*)
        { cout << "   'updateDisplacement()' is not defined for this Solid!\n\n"; return; }
        
        virtual void  resetMatrixAndVector()
        { cout << "   'resetMatrixAndVector()' is not defined for this Solid!\n\n"; return; }

        virtual void  AssembleElementVector(int ind, bool flag, double* rhs)
        { cout << "   'AssembleElementVector()' is not defined for this Solid!\n\n"; return; }
};






#endif
