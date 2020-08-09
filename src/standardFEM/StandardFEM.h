#ifndef incl_StandardFEM_CLASS_h
#define incl_StandardFEM_CLASS_h


#include "Solver.h"
#include "util.h"
#include "Domain.h"
#include "GeomDataLagrange.h"
#include "SolutionData.h"

#include "headersVTK.h"

#include "LagrangeElement.h"


class PropertyItem;


#define STANDARDFEM_ELEMENT_TYPE_NAMES {"LagrangeElem1Dbar2Node",\
                                    "EulerBernoulliBeamElement2D",     \
                                    "FrameElement2D",     \
                                    "ElementGeomExactTruss2D",     \
                                    "ElementGeomExactBeam2D", \
                                    "LagrangeElem2DPoissonTria3Node", \
                                    "LagrangeElem2DPoissonQuad4Node", \
                                    "LagrangeElem2DStructSolidTria3Node", \
                                    "LagrangeElem2DStructSolidQuad4Node", \
                                    "LagrangeElem2DStructSolidMixed", \
                                    "LagrangeElem2DBbarFbar", \
                                    "LagrangeElem2DStokesTria3Node", \
                                    "LagrangeElem2DStokesQuad4Node", \
                                    "LagrangeElem2DNavierStokesTria3Node", \
                                    "LagrangeElem2DNavierStokesQuad4Node", \
                                    "EulerBernoulliBeamElement3D", \
                                    "LagrangeElem3DPoissonTet4Node", \
                                    "LagrangeElem3DPoissonHex8Node", \
                                    "LagrangeElem3DStructSolidTet4Node", \
                                    "LagrangeElem3DStructSolidHex8Node", \
                                    "LagrangeElem3DStructSolidMixed", \
                                    "LagrangeElem3DBbarFbar", \
                                    "LagrangeElem3DStokesTet4Node", \
                                    "LagrangeElem3DStokesHex8Node", \
                                    "LagrangeElem3DNavierStokesTet4Node", \
                                    "LagrangeElem3DNavierStokesHex8Node", \
                                    "LagrangeElem2DStructSolidTria3NodeStab", \
                                    "LagrangeElem2DStructSolidQuad4NodeStab", \
                                    "LagrangeElem3DStructSolidTet4NodeStab", \
                                    "LagrangeElem3DStructSolidHex8NodeStab", \
                                    "MindlinPlateElement", \
                                    "KirchhoffPlateElement", \
                                    "LagrangeElem3DShellQuad4Node", \
                                    "ContactElement2D1nodedContactAlongXaxis", \
                                    "ContactElement2D1nodedContactAlongYaxis", \
                                    "ContactElement3D1nodedContactAlongXaxis", \
                                    "ContactElement3D1nodedContactAlongYaxis", \
                                    "ContactElement3D1nodedContactAlongZaxis",NULL}


class StandardFEM: public Domain
{
    //private:

    public:

        PhysicsType  PHYSICS_TYPE;
        SolverType   SOLVER_TYPE;

        int  DIM, ndof, lumpType, filecount, numProc, nElem_Constraint;
        int  nNode, nElem, npElem, IterNum, iterCount, totalDOF, fluidDOF, velDOF, presDOF;

        double  PENALTY_NUM, rNormPrev, rNorm,  ctimFactSolvUpdt, ctimCalcStiffRes, totalError, rhoInfty;

        bool  localStiffnessError, STAGGERED;
        char  VTKfilename[200];

        PetscInt  n_mpi_procs, this_mpi_proc;
        PetscInt  row_start, row_end, ndofs_local;
        PetscInt  elem_start, elem_end, nElem_local, nNode_local;

        PetscErrorCode ierr;

        vector<myPoint>  nodePosData;

        vector<vector<int> >  IEN, IEN2, ID, LM;
        vector<vector<int> >  boundaryNodes, forAssyMat;
        vector<vector<bool> >  NodeType;

        vector<int>  assy4r;

        VectorXd  soln, solnInit, reac;
        VectorXd  totalForce, totalMoment, centroid;

        vector<vector<double> >  DirichletBCs, NeumannBCs, DerivativeBCs, pointBCs, OutputData, nodeForcesData;

        vector<vector<int> >  elemConn, DDconn, DDconnE, DDconnEt, DDconnG, DDconnH, DDconnHt;

        vector<int>  node_map_new_to_old;
        vector<int>  node_map_old_to_new;

        vector<int>  dof_map_new_to_old;
        vector<int>  dof_map_old_to_new;

        LagrangeElement  **elems;

        MyString   anlySolnType;

        myPoint  geom, param, normal;

        SolutionData  SolnData;
        GeomDataLagrange  GeomData;

        SolverEigen *solverEigen;
        SolverPetsc *solverPetsc;

        vtkSmartPointer<vtkUnstructuredGrid>     uGridVTK;
        vtkSmartPointer<vtkPoints>               pointsVTK;
        vtkSmartPointer<vtkVertex>               vertexVTK;
        vtkSmartPointer<vtkLine>                 lineVTK;
        vtkSmartPointer<vtkQuad>                 quadVTK;
        vtkSmartPointer<vtkHexahedron>           hexVTK;

        vtkSmartPointer<vtkTriangle>             triaVTK;
        vtkSmartPointer<vtkPolygon>              polygonVTK;
        vtkSmartPointer<vtkTetra>                tetraVTK;
        vtkSmartPointer<vtkPyramid>              pyramidVTK;
        vtkSmartPointer<vtkWedge>                wedgeVTK;

        vtkSmartPointer<vtkIntArray>          procIdVTK;
        vtkSmartPointer<vtkFloatArray>          vecVTK, vecVTK2, scaVTK, scaVTK2, cellDataVTK, cellDataVTK2;
        vtkSmartPointer<vtkXMLUnstructuredGridWriter>  writerUGridVTK;

    public:

        StandardFEM();
        
        ~StandardFEM();

        ///////////////////////////////////////////////////////////
        //
        // DATA related member functions
        //
        ///////////////////////////////////////////////////////////

        void setDimension(int dd)
        {  ndm = DIM = dd;        }

        void setNdof(int dd)
        {  ndf = ndof = dd;        }

        void  setPhysicsTypetoSolid()
        {
          PHYSICS_TYPE = PHYSICS_TYPE_SOLID;
          SolnData.setPhysicsTypetoSolid();
          return;
        }

        void  setPhysicsTypetoFluid()
        {
          PHYSICS_TYPE = PHYSICS_TYPE_FLUID;
          SolnData.setPhysicsTypetoFluid();
          return;
        }

        void  setVTKfilename(char* ff)
        { std::strcpy(VTKfilename, ff);  return; }

        void  setNodes(const vector<myPoint>&  nodeVec)
        {
          nodePosData = nodeVec;
        }

        void  readNodes(ifstream& fname);

        void  readElementConnectivity(ifstream& fname);

        void  readElementProps(ifstream& fname);

        void  readMaterialProps(ifstream& fname);

        void  readInput(ifstream& fname);

        void  setElementConnectivity(const vector<vector<int> >&  nodeVec)
        {
          elemConn = nodeVec;
        }

        void  setNodeType(const vector<vector<bool> >&  nodeVec)
        {
          NodeType = nodeVec;
        }

        void  addElementProperties(const PropertyItem&  elmProp)
        {
          //ElemProp.add(elmProp);
          //ElemProp.add(new PropertyItem(ELEMENTTYPE));
          //ElemProp[ElemProp.n-1] = elmProp;
        }

        void  addMaterialProperties(const PropertyItem&  elmProp)
        {
          //MatlProp.add(elmProp);
          //MatlProp.add(new PropertyItem(MATERIAL));
          //MatlProp[MatlProp.n-1] = elmProp;
        }

        void  setDirichletBCs(const vector<vector<double> >&  nodeVec)
        {
          DirichletBCs = nodeVec;
        }

        void  setNodalForces(const vector<vector<double> >&  nodeVec)
        {
          nodeForcesData = nodeVec;
        }

        void  setControl(int tis1, double tol1, double rho1)
        {
          tis = tis1;
          rhoInfty = rho1;
          tol = tol1;

          SolnData.setTimeIncrementType(tis);
          SolnData.setSpectralRadius(rhoInfty);
        }

        virtual void printData(int, int);

        int  solveStep(int niter);

        ///////////////////////////////////////////////////////////
        //
        // PRE-PROCESSOR PHASE member functions
        //
        ///////////////////////////////////////////////////////////

        virtual void printComputerTime(bool reset = true, int detailFlg = 1);

        void assignBoundaryConditions();

        void  prepareElemProp();
        void  prepareMatlProp();

        virtual void prepareInteractions();

        virtual void prepareInputData();

        virtual void findMinMaxX(double*, double*, bool);

        int  findCellNumber(myPoint& geom);

        void  geometryToParametric(const myPoint& geom, myPoint& param);

        double  computeGeometry(const int dir, double param);

        void  computeGeometry(const myPoint& param, myPoint& geom);

        void  printInfo();

        virtual void readInputData(std::ifstream &, MyString &);

        void  plotGeom(int, bool, int, bool, int*);

        void  calcForceVector();

        //void  applyBoundaryConditions(int start, SparseMatrixXd& globalK, double* rhs);

        //void  applyBoundaryConditions();

        void  applyExternalForces();

        void  writeNodalData();

        void  writeReadResult(int, MyString &, int);

        ///////////////////////////////////////////////////////////
        //
        // SOLUTION PHASE member functions
        //
        ///////////////////////////////////////////////////////////


        virtual void setSolver(int, int *parm = NULL, bool cIO = false);

        int  prepareMatrixPattern();

        virtual int calcStiffnessAndResidual(int printRes=2, bool zeroMtx=true, bool zeroRes=true);

        virtual int factoriseSolveAndUpdate();

        virtual bool converged();

        virtual bool diverging(double);

        virtual void setTimeParam();

        virtual void timeUpdate();

        virtual void updateIterStep();

        virtual void reset();

        virtual void addExternalForces();

        void  applyBoundaryConditions();
    
        void  computeElementErrors(int);

        void  setInitialConditions();


        ///////////////////////////////////////////////////////////
        //
        // POST-PROCESSOR PHASE member functions
        //
        ///////////////////////////////////////////////////////////

        void  contourplot();

        void  postProcess(int, int, int, bool, double, double, int*);

        void  computeTotalBodyForce(int );

};







#include "DomainInlineFunctions.h"

define_reference_cast(standardFEM, StandardFEM)

define_isType(isStandardFEM, STANDARDFEM)


#endif














