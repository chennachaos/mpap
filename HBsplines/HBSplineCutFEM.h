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

// forAssyCutFEM   -  a vector<vector<int> > of size 'numDomains X gridDOF1*ndof' 
//                    the DOF that are not present in the current domain are given a value of '-1'
//                    This data is used for assembling element matrix and vector

// forAssyCutFEM2  -  a vector<vector<int> > of size 'numDomains X domainTotalDOF[]'. 
//                    This vector stores the global positions of each active DOF for 
//                    the current domain.
//                    This data is used for extracting solution for each domain from
//                    the global solution vector


class HBSplineCutFEM: public HBSplineBase
{
    private:

        int  numDomains, fluidDOF, solidDOF, CUTCELL_INTEGRATION_TYPE;

        double  GHOST_PENALTY;

        vector<int>  domainTotalDOF, domainStartDOF, cutCellIds, fluidElementIds;

        vector<bool>  domainInclYesNo;

        vector<int>  forAssyCutFEM, forAssyCutFEMprev, forAssyCutFEM2;

        vector<double>  cutFEMparams;

        VectorXd  totalForce, slnTemp,  slnTempPrev, slnTempPrev2, slnTempCur;

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

};







#include "DomainInlineFunctions.h"

define_reference_cast(hbsplineCutFEM, HBSplineCutFEM)

define_isType(isHBSplineCutFEM, HBSPLINECUTFEM)



#endif














