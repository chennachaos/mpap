#ifndef incl_HBSplineCutFEM_CLASS_h
#define incl_HBSplineCutFEM_CLASS_h


#include "headersBasic.h"
#include "HBSplineBase.h"


////////////////////////////////

// numDomains      -  total number of domains in the current problem

// domainTotalDOF  -  a vector<int> of size 'numDomains' containing total number of DOF in each domain

// domainStartDOF  -  a vector<int> of size 'numDomains' containing starting number for 
//                    element matrix assembly for each domain

// a domain can be active or inactive.
// For an 'active' domain all the compuatations are to be carried out
// For an 'inactive' domain no compuatations are to be carried out
// An 'inactive' domain represents holes in the global domain

// domainInclYesNo -  a vector<bool> of size 'numDomains' 
//                    0 -> the domain in inactive
//                    1 -> the domain in active

// grid_to_cutfem_BF    -  a vector<vector<int> > of size 'numDomains X gridDOF1' 
//                    the DOF that are not present in the current domain are given a value of '-1'
//                    This data is used for assembling element matrix and vector.
//                    For fluid flow this is just a vector as there is only one active domain.

// grid_to_cutfem_DOF    -  a vector<vector<int> > of size 'numDomains X gridDOF1*ndof' 
//                    the DOF that are not present in the current domain are given a value of '-1'
//                    This data is used for assembling element matrix and vector.
//                    For fluid flow this is just a vector as there is only one active domain.

// cutfem_to_grid_BF  -    a vector<vector<int> > of size 'numDomains X domainTotalDOF[]'. 
//                        This vector stores the global positions of each active DOF for 
//                        the current domain.
//                        This data is used for extracting solution for each domain from
//                        the global solution vector
//                        For fluid flow this is just a vector as there is only one active domain.

// cutfem_to_grid_DOF  -   a vector<vector<int> > of size 'numDomains X domainTotalDOF[]'. 
//                        This vector stores the global positions of each active DOF for 
//                        the current domain.
//                        This data is used for extracting solution for each domain from
//                        the global solution vector
//                        For fluid flow this is just a vector as there is only one active domain.

// grid_to_proc_BF[ii]     - node_map_old_to_new[grid_to_cutfem_BF[ii]]
// grid_to_proc_DOF[ii]    - dof_map_old_to_new[grid_to_cutfem_DOF[ii]]

// proc_to_grid_BF[ii]     - cutfem_to_grid_BF[node_map_new_to_old[ii]]
// proc_to_grid_DOF[ii]    - cutfem_to_grid_DOF[dof_map_new_to_old[ii]]


class HBSplineCutFEM: public HBSplineBase
{
    private:

        int  numDomains, fluidDOF, solidDOF, CUTCELL_INTEGRATION_TYPE;

        vector<int>  domainTotalDOF, domainStartDOF, cutCellIds, fluidElementIds;

        vector<bool>  domainInclYesNo;

        vector<int>  grid_to_cutfem_BF, grid_to_cutfem_BFprev, grid_to_cutfem_DOF, grid_to_cutfem_DOFprev;
        vector<int>  cutfem_to_grid_BF, cutfem_to_grid_DOF;
        vector<int>  grid_to_proc_BF, grid_to_proc_DOF, proc_to_grid_BF, proc_to_grid_DOF;

        vector<double>  cutFEMparams;

        VectorXd  slnTemp,  slnTempPrev, slnTempPrev2, slnTempCur;

        vector<set<int> >  DDconnLoc;

    public:

        HBSplineCutFEM();

        virtual ~HBSplineCutFEM();

        virtual  void  prepareInputData();

        virtual  void  printInfo();

        virtual  void  readInputData(std::ifstream &, MyString &);

        virtual  void  writeReadResult(int, MyString &);

        virtual  void  plotGeom(int, bool, int, bool, int*);

        void  plotGeomSubTrias1D(int, bool, int, bool, int*);
        void  plotGeomSubTrias2D(int, bool, int, bool, int*);
        void  plotGeomSubTrias3D(int, bool, int, bool, int*);

        void  plotGeomAdapIntegration1D(int, bool, int, bool, int*);
        void  plotGeomAdapIntegration2D(int, bool, int, bool, int*);
        void  plotGeomAdapIntegration3D(int, bool, int, bool, int*);

        virtual  void  plotGaussPointsElement();

        virtual  void  plotGaussPointsDirichletBoundary();

        virtual  void  plotGaussPointsNeumannBoundary();

        virtual  void  setTimeParam();

        virtual  void  timeUpdate();

        virtual  void  updateIterStep();

        virtual  void  reset();

        virtual  int  prepareMatrixPattern();

        virtual  void  prepareMatrixPatternPostProcess();

        virtual  int  calcStiffnessAndResidual(int printRes=2, bool zeroMtx=true, bool zeroRes=true);

        void  applyInterfaceTerms2D();
        void  applyInterfaceTerms3D();

        void  applyGhostPenalty2D();
        void  applyGhostPenalty3D();

        int  setCoveringUncovering();

        int  setCoveringUncovering1D();
        int  setCoveringUncovering2D();
        int  setCoveringUncovering3D();

        virtual  void  solveSolidProblem();

        virtual  int  factoriseSolveAndUpdate();

        virtual  void  addExternalForces();

        virtual  void  applyBoundaryConditions();

        virtual  void  computeElementErrors(int);

        virtual  void  setInitialConditions();

        virtual  void  computeConditionNumber();

        virtual  void  computeTotalForce(int index);

        void  computeTotalForce2D(int index);

        void  computeTotalForce3D(int index);

        virtual  void  postProcessFlow(int, int, int, bool, double, double, int*);

        void  postProcessSubTrias1D(int, int, int, bool, double, double, int*);
        void  postProcessSubTrias2D(int, int, int, bool, double, double, int*);
        void  postProcessSubTrias3D(int, int, int, bool, double, double, int*);

        void  postProcessAdapIntegration1D(int, int, int, bool, double, double, int*);
        void  postProcessAdapIntegration2D(int, int, int, bool, double, double, int*);
        void  postProcessAdapIntegration3D(int, int, int, bool, double, double, int*);

        virtual  void  computeTotalBodyForce(int );

        void  prepareCutElements();

        void  prepareCutElementsSubTrias1D();
        void  prepareCutElementsSubTrias2D();
        void  prepareCutElementsSubTrias3D();

        void  prepareCutElementsAdapIntegration1D();
        void  prepareCutElementsAdapIntegration2D();
        void  prepareCutElementsAdapIntegration3D();

        //void  triaVTK(int kk, vtkIdType pt);

        void  computeGaussPoints();

        void  computeGaussPointsSubTrias1D();
        void  computeGaussPointsSubTrias2D();
        void  computeGaussPointsSubTrias3D();

        void  computeGaussPointsAdapIntegration1D();
        void  computeGaussPointsAdapIntegration2D();
        void  computeGaussPointsAdapIntegration3D();

        //virtual  int  solveFluidProblem();

        virtual  int  deallocatePetscObjects();

};







#include "DomainInlineFunctions.h"

define_reference_cast(hbsplineCutFEM, HBSplineCutFEM)

define_isType(isHBSplineCutFEM, HBSPLINECUTFEM)



#endif














