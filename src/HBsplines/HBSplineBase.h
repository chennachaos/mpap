
#ifndef incl_HBSplineBase_CLASS_h
#define incl_HBSplineBase_CLASS_h


#include "Solver.h"
#include "TreeNode.h"
#include "Domain.h"
#include "Definitions.h"
#include "DomainTree.h"
#include "headersVTK.h"
#include "headersEigen.h"
#include "ImmersedSolid.h"
#include "SolverEigen.h"
#include "SolverPetsc.h"


using namespace std;


class ContactElementPointToPoint2D;
class SolutionData;
class GeomDataHBSplines;



#define HBS_ELEMENT_TYPE_NAMES {"HBSplineElem1DAdvectionDiffusion",  \
                                "HBSplineElem2DAdvectionDiffusion", \
                                "HBSplineElem2DStokes", \
                                "HBSplineElem2DNavierStokes3dof", \
                                "HBSplineElem2DNavierStokes4dof", \
                                "HBSplineElem2DHeatTransfer", \
                                "HBSplineElem2DTempCoupled4dof", NULL}


class HBSplineBase: public Domain
{
    //typedef TreeNode<2>  node;
    //typedef  point3d  myPoint;

    public:

        int  DIM, CURRENT_LEVEL, MAX_LEVEL, numActiveElem, ndof, numActiveBasis, filecount, LSFEM_FLAG, numProc;
        int  degree[3], nlbf[3], nelem[3], totnumel, IterNum, iterCount;
        int  gridBF1, gridBF2, gridBFtotal, totalDOF;
        int  nImmSolids, SOLVER_TYPE, slv_type, SOLID_SOLVED;

        PetscInt  nElem, nNode, nElem_local, nNode_local;
        PetscInt  n_mpi_procs, this_mpi_proc;
        PetscInt  bfs_start, bfs_end, bfs_local;
        PetscInt  row_start, row_end, ndofs_local;
        PetscInt  elem_start, elem_end;

        PetscInt  *colTemp;
        PetscScalar  *arrayTemp;

        PetscErrorCode ierr;

        double  origin[3], gridLEN[3], knotsLeft[3], knotsRight[3];
        double  rNormPrev, rNorm,  ctimFactSolvUpdt, ctimCalcStiffRes, totalError, rhoInfty;

        vector<myPoint>  gridVertices;

        bool  STATUS_BCS, localStiffnessError, PERIODIC_BCS, CREATE_POSTPROCESS_GRID;
        bool  GRID_CHANGED, IB_MOVED, STAGGERED, STABILISED;

        vector<int>  VacantBFs, nodes2divide, elemsToRefine, elemsToUnRefine, activeElements;

        vector<double>  refinementData, fluidProps, stagParams;
        vector<vector<int> >  boundaryNodes, contElemData;
        vector<vector<int> >  NodeNumsAtLevel;
        vector<vector<int> >  DDconn;

        vector<vector<double> >  DirichletBCs, NeumannBCs, DerivativeBCs, pointBCs, Iconds;
        vector<vector<double> >  refineLimitVals, FluidOutputData;

        vector<node*>  elems;
        //typedef boost::ptr_vector<base> container;
        //boost::ptr_vector<node>  elems;
        //std::vector<std::unique_ptr<node> > elems;
        //vector< std::tr1::shared_ptr<node> > elems;

        vector<int>  node_map_new_to_old, node_map_old_to_new;
        vector<int>  dof_map_new_to_old, dof_map_old_to_new;

        MyString   anlySolnType;

        List<PropertyItem>  ElemProp, MatlProp;

        node  *root;

        myPoint  geom, param, normal;

        vector<ImmersedSolid*>  ImmersedBodyObjects;
        vector<ContactElementPointToPoint2D*>  contactElementObjects;

        SparseMatrixXd  globalK2;

        VectorXd   soln, solnInit, totalForce;

        SolutionData  SolnData;

        GeomDataHBSplines  GeomData;

        SolverEigen  *solverEigen;

        SolverPetsc  *solverPetsc;

        BiCGSTAB<SparseMatrixXd, IncompleteLUT<double> > solver3;

        vtkSmartPointer<vtkUnstructuredGrid>     uGridVTK;
        vtkSmartPointer<vtkPoints>               pointsVTK;
        vtkSmartPointer<vtkVertex>               vertexVTK;
        vtkSmartPointer<vtkLine>                 lineVTK;
        vtkSmartPointer<vtkQuad>                 quadVTK;
        vtkSmartPointer<vtkTetra>                tetVTK;
        vtkSmartPointer<vtkHexahedron>           hexVTK;

        vtkSmartPointer<vtkFloatArray>          vecVTK, vecVTK2, scaVTK, scaVTK2, cellDataVTK, cellDataVTK2;
        vtkSmartPointer<vtkXMLUnstructuredGridWriter>  writerUGridVTK;

        vtkSmartPointer<vtkMergePoints>   mergePoints;
        vtkSmartPointer<vtkPoints>               pointsVTKfluidgrid;
        vtkSmartPointer<vtkPolyData>            pointsPolydata ;
        vtkSmartPointer<vtkUnstructuredGrid>     uGridVTKfluid;

    public:

        HBSplineBase();

        virtual ~HBSplineBase();

        ///////////////////////////////////////////////////////////
        //
        // DATA related member functions
        //
        ///////////////////////////////////////////////////////////

        int  getMaxLevel()
        {  return MAX_LEVEL;  }

        bool isLeaf(int nodenum)
        {  return false;      }

        bool isActive(int nodenum)
        {  return isLeaf(nodenum);   }
        
        int getDimension()
        {  return DIM;        }

        int* getDegree()
        {  return  degree;        }

        int getNumberOfActiveElements()
        {  return  numActiveElem;        }

        int getNumberOfActiveBasisFunctions()
        {  return numActiveBasis;        }

        node*  getRoot()
        {  return  root;        }

        void setDimension(int dd)
        {  ndm = DIM = dd;        }

        void setNdof(int dd)
        {  ndf = ndof = dd;        }

        void  setOrigin(double x0=0.0, double y0=0.0, double z0=0.0)
        {
          origin[0] = x0;
          origin[1] = y0;
          origin[2] = z0;
        }

        void  setGridDimensions(double xl=1.0, double yl=1.0, double zl=1.0)
        {
          gridLEN[0] = xl;
          gridLEN[1] = yl;
          gridLEN[2] = zl;
        }

        void  setBSplineDegree(int p1=1, int p2=0, int p3=0)
        {
          degree[0] = p1;
          degree[1] = p2;
          degree[2] = p3;
        }

        void  setNumberOfElements(int n1=1, int n2=1, int n3=1)
        {
          nelem[0] = n1;
          nelem[1] = n2;
          nelem[2] = n3;
        }

        void  addDirichletBCs(int side, int dof, double val, double PEN, int type, double fact)
        {
          int aa = DirichletBCs.size();
          vector<double>  vectemp(6);
          DirichletBCs.push_back(vectemp);

          DirichletBCs[aa][0] = side;
          DirichletBCs[aa][1] = dof;
          DirichletBCs[aa][2] = val;
          DirichletBCs[aa][3] = PEN;
          DirichletBCs[aa][4] = type;
          DirichletBCs[aa][5] = fact;
        }

        void  addNeumannBCs(int side, int dof, double val, double PEN, int type, double fact)
        {
          int aa = NeumannBCs.size();
          vector<double>  vectemp(6);
          NeumannBCs.push_back(vectemp);

          NeumannBCs[aa][0] = side;
          NeumannBCs[aa][1] = dof;
          NeumannBCs[aa][2] = val;
          NeumannBCs[aa][3] = PEN;
          NeumannBCs[aa][4] = type;
          NeumannBCs[aa][5] = fact;
        }

        void  addPointBCs(double xx, double yy, double zz, int dof, double val, double PEN)
        {
          int aa = pointBCs.size();
          vector<double>  vectemp(6);
          pointBCs.push_back(vectemp);

          pointBCs[aa][0] = xx;
          pointBCs[aa][1] = yy;
          pointBCs[aa][2] = zz;
          pointBCs[aa][3] = dof;
          pointBCs[aa][4] = val;
          pointBCs[aa][5] = PEN;
        }

        void setFluidProperties(const vector<double>&  vectemp)
        {
          fluidProps = vectemp;
        }

        void setControl(int tis1, double tol1, double rho1)
        {
           tis = tis1;
           td[0] = rho1;
           tol = tol1;
        }

        void  subDivide(int nodenum);

        void  addNeighbourElements(int kk);
        void  addNeighbourElements1D(int kk);
        void  addNeighbourElements2D(int kk);
        void  addNeighbourElements3D(int kk);
        
        void  unRefine(int nodenum);

        ///////////////////////////////////////////////////////////
        //
        // PRE-PROCESSOR PHASE member functions
        //
        ///////////////////////////////////////////////////////////

        void assignBoundaryConditions();

        void processBoundaryConditionsRefinedLevels();

        virtual void prepareInteractions();

        virtual void prepareInputData();

        virtual void findMinMaxX(double*, double*, bool);

        int  findCellNumber(const myPoint& geom);

        void  geometryToParametric(const myPoint& geom, myPoint& param);

        double  computeGeometry(const int dir, double param);

        void  computeGeometry(const myPoint& param, myPoint& geom);

        void  buildBase();
        void  buildBase1D();
        void  buildBase2D();
        void  buildBase3D();

        virtual void  printInfo();

        void  printSelf()
        {
          for(int ii=0;ii<elems.size();ii++)
            elems[ii]->printSelf();
        }

        void  printNodeInfo(int ind)
        {
            assert( (ind) < elems.size() );
            elems[ind]->printSelf();
        }
        
        void   refine(int kk);
        
        void   applyRefinementProcess();
        
        void   addGhostNodes(node* nd, int direction);

        void   algorithm1(int lev);
        void   algorithm2(int lev);
        void   algorithm3(int lev);

        void   algorithm1_1D(int lev);
        void   algorithm2_1D(int lev);
        void   algorithm3_1D(int lev);

        void   algorithm1_2D(int lev);
        void   algorithm2_2D(int lev);
        void   algorithm3_2D(int lev);

        void   algorithm1_3D(int lev);
        void   algorithm2_3D(int lev);
        void   algorithm3_3D(int lev);

        virtual void readInputData(std::ifstream &, MyString &)
        { cout << " readInputData() ... is not defined for the class... " << endl; return; }

        virtual void  plotGeom(int, bool, int, bool, int*)
        { cout << " plotGeom() ... is not defined for the class... " << endl; return; }

        virtual void  plotGaussPoints()
        { cout << " plotGaussPoints() ... is not defined for the class... " << endl; return; }

        void  immersedBoundaryConditions2D1();
        void  immersedBoundaryConditions2D2();
        void  immersedBoundaryConditions2D3();

        void  createImmersedBoundaryPoints();

        virtual  void  solveSolidProblem()
        { cout << " solveSolidProblemwriteReadResult() ... is not defined for the class... " << endl; return; }

        void  updateImmersedPointPositions();

        virtual void writeNodalData();

        void  writeImmersedSolidOutput();
        void  writeFluidOutput();

        void  refinementforHemkerProblem();
        void  refinementforAdvDiff1D();
        void  refinementforAdvDiff2D();

        void  pointBasedRefinement(int kk);
        void  limitBasedRefinement(int kk);

        virtual  void  writeReadResult(int, MyString &)
        { cout << " writeReadResult() ... is not defined for the class... " << endl; return; }

        ///////////////////////////////////////////////////////////
        //
        // SOLUTION PHASE member functions
        //
        ///////////////////////////////////////////////////////////

        virtual void printData(int, int);

        virtual void setSolver(int, int *parm = NULL, bool cIO = false);

        virtual bool converged();

        virtual bool diverging(double);

        virtual void setTimeParam()
        { cout << " setTimeParam() ... is not defined for the class... " << endl; return; }

        virtual void timeUpdate()
        { cout << " timeUpdate() ... is not defined for the class... " << endl; return; }

        virtual void updateIterStep()
        { cout << " updateIterStep() ... is not defined for the class... " << endl; return; }

        virtual void reset()
        { cout << " reset() ... is not defined for the class... " << endl; return; }

        virtual void printComputerTime(bool reset = true, int detailFlg = 1);

        virtual  int  prepareMatrixPattern()
        { cout << " prepareMatrixPattern() ... is not defined for the class... " << endl; return 1; }

        virtual void prepareMatrixPatternPostProcess()
        { cout << " prepareMatrixPatternPostProcess() ... is not defined for the class... " << endl; return; }

        virtual int calcStiffnessAndResidual(int printRes=2, bool zeroMtx=true, bool zeroRes=true)
        { cout << " calcStiffnessAndResidual() ... is not defined for the class... " << endl; return -1; }

        virtual int factoriseSolveAndUpdate()
        { cout << " factoriseSolveAndUpdate() ... is not defined for the class... " << endl; return -1; }

        virtual void addExternalForces()
        { cout << " addExternalForces() ... is not defined for the class... " << endl; return; }

        virtual void  applyBoundaryConditions()
        { cout << " applyBoundaryConditions() ... is not defined for the class... " << endl; return; }

        virtual void  computeElementErrors(int)
        { cout << " computeElementErrors() ... is not defined for the class... " << endl; return; }

        virtual void  performAdaptiveRefinement(double eTol);

        virtual void  setInitialConditions();

        virtual void  computeConditionNumber();

        ///////////////////////////////////////////////////////////
        //
        // POST-PROCESSOR PHASE member functions
        //
        ///////////////////////////////////////////////////////////

        virtual  void  computeTotalForce(int index)
        { cout << " computeTotalForce() ... is not defined for the class... " << endl; return; }

        void  printResultAtPoint(int, double, double, double);

        virtual void  postProcessFlow(int, int, int, bool, double, double, int*)
        { cout << " postProcessFlow() ... is not defined for the class... " << endl; return; }

        void  computeTotalBodyForce(int );

        virtual  int  solveFluidProblem();

        virtual  int  fsi_staggered_force_predictor(int max_iter, double tol_local);

        virtual  int  fsi_staggered_displacement_predictor(int max_iter, double tol_local);

        virtual  int  fsi_monolithic_fixedpoint_forcePred(int max_iter, double tol_local);

        virtual  int  fsi_monolithic_fixedpoint_dispPred(int max_iter, double tol_local);

        //virtual  int  deallocatePetscObjects();

};





#endif









