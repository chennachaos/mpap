#ifndef incl_ImmersedSolid_h
#define incl_ImmersedSolid_h

#include "SolutionData.h"
#include "GeomDataLagrange.h"
#include "AABB.h"
#include "Ray.h"
#include "myPoly.h"
#include "headersVTK.h"
#include  "SolverEigen.h"
#include  "SolverPetsc.h"

#include "TimeFunctionCore.h"
#include "myCGALroutines.h"

/*
#include <list>
//#include <cfloat>
//#define DBL_MAX    1.7976931348623157E+308

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

//typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::FT CGAL_FT;
typedef K::Ray_3 CGAL_Ray;
typedef K::Line_3 CGAL_Line;
typedef K::Point_3 CGAL_Point;
typedef K::Segment_3 CGAL_Segment;
typedef K::Triangle_3 CGAL_Triangle;
typedef std::list<CGAL_Triangle>::iterator CGAL_Iterator;
typedef CGAL::AABB_triangle_primitive<K, CGAL_Iterator> CGAL_Primitive;
typedef CGAL::AABB_traits<K, CGAL_Primitive> CGAL_AABB_triangle_traits;
typedef CGAL::AABB_tree<CGAL_AABB_triangle_traits> CGAL_Tree;
*/

using namespace myGeom;

class  ImmersedIntegrationElement;


enum  {BC_ENFORCE_TYPE_LAGRANGE=0, BC_ENFORCE_TYPE_PENALTY};


class ImmersedSolid
{
    public:

        static  int  count;

        int  id, DIM, type, localStiffnessError, filecount, iterCount, tis, ndof, solverOK ;
        int  nElem, nImmInt, totalDOF, npElem, nNode, nsize, nElem_Constraint;
        int  matId, elemId, BC_ENFORCE_TYPE;
        int  n_mpi_procs, this_mpi_proc;

        double  rNorm, rNormPrev, tol, dt, rho, beta, ctimCalcStiffRes, ctimFactSolvUpdt, PENALTY, NitscheFact;

        bool firstIter, STAGGERED, isNitsche, PRESC_MOTION;

        vector<int>  assy4r;
        vector<double>  preLoad, initForcePred;
        vector<vector<int> >  forAssyMat, forAssyCoupledHorz, forAssyCoupledVert;
        vector<vector<int> >  OutputData;
        vector<vector<double> > rigidBodyMotionLimits;
        vector<TimeFunctionCore> PrescMotionTimeFuncs;

        VectorXd  soln, totalForce, fluidAcce, fluidAccePrev, fluidAcceCur;
        myPoint  centroid, pivotpoint;

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

        // CGAL related
        vector<CGAL_Point>   pointsCGAL;

        list<CGAL_Triangle>  trianglesCGAL;
        CGAL_Tree treeCGAL;

        vector<int>  facesCGAL;
        CGAL_Polyhedron  polyhedronTemp;
        CGAL_Point_inside  *point_inside_tester;

        //////////////////////////////////////////////////////
        // member functions
        //////////////////////////////////////////////////////

        ImmersedSolid();

        virtual ~ImmersedSolid();

        int getID()
        {  return id; }

        void  setDimension(int dd)
        { DIM = dd; }

        void  setNumberOfNodes(int nn)
        {  nNode = nn; }

        int  getNumberOfNodes()
        { return nNode; } 

        void  setNumberOfElements(int nn)
        {  nElem = nn; }

        int  getNumberOfElements()
        {  return nElem; }

        void  setNumberOfNodesPerElement(int nn)
        {  npElem = nn; }

        int getTotalDOF()
        { return  totalDOF; }

        void  setTolerance(double tt)
        { tol = tt; }

        void  setPenaltyParameter(double tt)
        { PENALTY = tt; }

        double  getPenaltyParameter()
        { return PENALTY; }

        void  setNitscheFact(double tt)
        { NitscheFact = tt; }

        double  getNitscheFact()
        { return NitscheFact; }

        void  setNitscheFlag(bool tt)
        { isNitsche = tt; }

        bool  getNitscheFlag()
        { return isNitsche; }

        void setTimeIncrementType(int ttt)
        {  SolnData.setTimeIncrementType(ttt); }

        void setSpectralRadius(double ttt)
        {  SolnData.setSpectralRadius(ttt); }

        void setInitDisplacement(int ind, int dir, double ttt)
        {  SolnData.var1[ind*ndof+dir] = ttt; }

        void setInitVelocity(int ind, int dir, double ttt)
        {  SolnData.var1Dot[ind*ndof+dir] = ttt; }

        void setInitAcceleration(int ind, int dir, double ttt)
        {  SolnData.var1DotDot[ind*ndof+dir] = ttt; }

        double getDisplacement(int ind, int dir)
        {  return  SolnData.var1[ind*ndof+dir]; }

        double getDisplacementCur(int ind, int dir)
        {  return  SolnData.var1Cur[ind*ndof+dir]; }

        double getVelocity(int ind, int dir)
        {  return  SolnData.var1Dot[ind*ndof+dir]; }

        double getVelocityCur(int ind, int dir)
        {  return  SolnData.var1DotCur[ind*ndof+dir]; }

        double getAcceleration(int ind, int dir)
        {  return  SolnData.var1DotDot[ind*ndof+dir]; }

        double getAccelerationCur(int ind, int dir)
        {  return  SolnData.var1DotDotCur[ind*ndof+dir]; }

        double getForce(int ind, int dir)
        {  return  SolnData.force[ind*ndof+dir]; }

        double getForceCur(int ind, int dir)
        {  return  SolnData.forceCur[ind*ndof+dir]; }

        void  setBoundaryConditionType(int tt)
        {
          if(tt)
            BC_ENFORCE_TYPE = BC_ENFORCE_TYPE_LAGRANGE;
          else
            BC_ENFORCE_TYPE = BC_ENFORCE_TYPE_PENALTY;
        }

        bool isBoundaryConditionTypeLagrange()
        {
          return (BC_ENFORCE_TYPE == BC_ENFORCE_TYPE_LAGRANGE);
        }

        void  setMaterialID(int tt)
        { matId = tt; }

        int  getMaterialID()
        { return matId; }

        void  setElementID(int tt)
        { elemId = tt; }

        int  getElementID()
        { return elemId; }

        virtual bool isRigidBody()
        { cout << "   'isRigidBody()' is not defined for this Solid!\n\n"; return true; }

        virtual bool isFlexibleBody()
        { cout << "   'isFlexibleBody()' is not defined for this Solid!\n\n"; return true; }

        virtual  void  setNodalPositions(vector<vector<double> >&  vectemp) = 0;

        virtual void  setSolidElements(vector<vector<int> >& datatemp)
        { cout << "   'SetSolidElements' is not defined for this Solid!\n\n"; return; }

        virtual void  setImmersedElemActiveFlag(vector<int>& datatemp);

        virtual void  setImmersedIntegrationElements(vector<vector<int> >& datatemp);

        virtual  void  setImmersedFaces();

        virtual  void  updateImmersedFaces();

        void  adjustBoundaryPoints(double* minVal, double* maxVal);

        virtual void  setDataForOutput(vector<vector<int> >& vectemp);

        virtual void  setBoundaryConditions(vector<vector<double> >& vectemp)
        { cout << "   'SetBoundaryConditions' is not defined for this Solid!\n\n"; return; }

        virtual void  printSelf()
        { cout << "   'printSelf()' is not defined for this Solid!\n\n"; return; }

        virtual void  initialise()
        { cout << "   'initialise()' is not defined for this Solid!\n\n"; return; }

        virtual void  setMass(vector<double>&)
        { cout << "   'setMass()' is not defined for this Solid!\n\n"; return; }

        virtual void  setDamping(vector<double>&)
        { cout << "   'setDamping()' is not defined for this Solid!\n\n"; return; }

        virtual void  setStiffness(vector<double>&)
        { cout << "   'setStiffness()' is not defined for this Solid!\n\n"; return; }

        virtual void  setBoundaryConditions(vector<int>& vectemp)
        { cout << "   'SetBoundaryConditions' is not defined for this Solid!\n\n"; return; }

        virtual void  setPrescribedMotion(vector<vector<double> >& vectemp)
        { cout << "   'SetPrescribedMotion' is not defined for this Solid!\n\n"; return; }

        virtual void  setPivotpoint(vector<double>&)
        { cout << "   'setPivotpoint' is not defined for this Solid!\n\n"; return; }

        virtual void  setPreload(vector<double>&)
        { cout << "   'setPreload' is not defined for this Solid!\n\n"; return; }

        virtual void  setInitialForcePredictor(vector<double>&)
        { cout << "   'setInitialForcePredictor' is not defined for this Solid!\n\n"; return; }

        virtual void  setRigidBodyMotionLimits(vector<vector<double> >&)
        { cout << "   'setRigidBodyMotionLimits' is not defined for this Solid!\n\n"; return; }

        virtual void  setSolver(int slv = 1, int *parm = NULL, bool cIO = false)
        { cout << "   'setSolver()' is not defined for this Solid!\n\n"; return; }

        virtual void  computeInitialAcceleration()
        { cout << "   'computeInitialAcceleration()' is not defined for this Solid!\n\n"; return; }

        virtual int  calcStiffnessAndResidual(int solver_type=1, bool zeroMtx=true, bool zeroRes=true)
        { cout << "   'calcStiffnessAndResidual()' is not defined for this Solid!\n\n"; return 0; }

        virtual int  applyBoundaryConditions(int start1, int start2, SparseMatrixXd& globalK, double* rhs)
        { cout << "   'applyBoundaryConditions()' is not defined for this Solid!\n\n"; return -1; }

        virtual int  applyExternalForces()
        { cout << "   'applyExternalForces()' is not defined for this Solid!\n\n"; return -1; }

        virtual int  factoriseSolveAndUpdate()
        { cout << "   'factoriseSolveAndUpdate()' is not defined for this Solid!\n\n"; return 0; }

        virtual void  setTimeParam();

        virtual bool  diverging(double factor);

        virtual bool  converged();

        virtual void  solveTimeStep()
        { cout << "   'solveTimeStep()' is not defined for this Solid!\n\n"; return; }

        virtual void  timeUpdate();

        virtual void  updateIterStep();

        virtual void  reset();

        void  computeCentroid(int index);

        myPoint&  getCentroid(int index)
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

        virtual int assembleGlobalMatrixAndVector(int ind1, int ind2, SparseMatrixXd& mtx, double* rhs)
        { cout << "   'assembleGlobalMatrixAndVector()' is not defined for this Solid!\n\n"; return 0; }

        virtual int assembleGlobalMatrixAndVectorCutFEM(int ind1, int ind2, SolverPetsc* solverTemp)
        { cout << "   'assembleGlobalMatrixAndVectorCutFEM()' is not defined for this Solid!\n\n"; return 0; }

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

        virtual void  assembleElementVector(int ind, bool flag, double* rhs)
        { cout << "   'assembleElementVector()' is not defined for this Solid!\n\n"; return; }

        virtual  void  writeResult(ofstream&)
        { cout << " writeResult() ... is not defined for the class... " << endl; return; }

        virtual  void  readResult(ifstream&)
        { cout << " readResult() ... is not defined for the class... " << endl; return; }

        virtual  void  perform_Aitken_accelerator_force()
        {   SolnData.perform_Aitken_accelerator_force(); return; }

        virtual  void  perform_Aitken_accelerator_displacement()
        {   SolnData.perform_Aitken_accelerator_displacement(); return; }

        void  predict_force();
        void  interpolate_force();

        void  predict_displacement();
        void  interpolate_displacement();

        void  initialise_solid_state();
};






#endif
