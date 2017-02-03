
#ifndef incl_TreeNode_h
#define incl_TreeNode_h

#include "headersEigen.h"

#include "AdaptiveOctree.h"
#include "AdaptiveBinarytree.h"
#include "GeomDataHBSplines.h"
#include "SolutionData.h"
#include "AABB.h"

#include <vector>

class  myDataIntegrateCutFEM;
class GaussQuadrature;


using namespace myGeom;
using std::vector;


template <int DIM>
class TreeNode
{
    typedef TreeNode<DIM>*  TreeNode_PTR;

    //typedef Matrix<double, DIM, 1>  eigenPoint_type;

    private:

        static int  nodecount, MAX_LEVEL;

        int  level, id, subdomId, NUM_CHILDREN, NUM_NEIGHBOURS, nlbf[DIM], nlbf2[DIM], degree[DIM];

        int  totnlbf, totnlbf2, nsize, nsize2, ndof, counter;
        NodeOrientation  orientation;
        //elemTypeCutFEM  elemType;

        double  *elmDat, *matDat, JacMultElem, elemError, volume;
        double  knots[3][4];
        
        AABB  bbox;
        
        //eigenPoint_type  knotBegin, knotEnd;
        myPoint  knotBegin, knotEnd, knotIncr;

        bool FuncFlag, GHOST_FLAG, ACTIVE_FLAG, PROCESSED;
    
        TreeNode_PTR  *child, parent, *neighbours;

        //int domainNum;

    public:

        vector<int>  domNums, QuadratureDomNums;

        vector<vector<double> > DirichletData, NeumannData, DerivativeBCData;

        vector<int>  GlobalBasisFuncs, LocalBasisFuncs, LocalBasisFuncsPrev, forAssyVec, forAssyVec2;

        SolutionData  *SolnData;
        GeomDataHBSplines  *GeomData;

        MatrixXd  SubDivMat, Klocal;

        VectorXd  Flocal;
        
        GaussQuadrature  Quadrature;
        vector<GaussQuadrature>  BoundaryQuadrature;

        vector<myPoly*>  subTrias;

        //typedef  AdaptiveOctree<DIM>*  adapTreePTR;

        typedef  AdaptiveBinarytree<DIM>*  adapTreePTR;

        //AdaptiveOctree<DIM>  *adapIntegNode1;
        //AdaptiveBinarytree<DIM>  *adapIntegNode2;
        adapTreePTR  adapIntegNode;

        TreeNode();

        TreeNode(int lev=0);

        ~TreeNode();

        int  GetLevel()
        {  return level; }

        void SetLevel(int lev)
        {  level = lev; }

        bool IsLeaf()
        {  return (child == NULL); }

        bool IsActive()
        {  return  ACTIVE_FLAG; }

        bool IsProcessed()
        {  return PROCESSED;}

        int GetID()
        {  return id;}

        int GetDimension()
        {  return DIM;}

        void  Setndof(int nd)
        { ndof = nd; }

        int  Getndof()
        { return ndof; }

        int* GetDegree()
        {  return  degree; }

        double  GetError()
        {  return  elemError; }

        double  GetJacMultElement()
	{ return JacMultElem; }

        double  GetVolume();
        
        double  GetVolumeGaussPoints(int domTemp);

        int GetNsize()
        { return nsize2; }

        void  set_subdomain_id(int sid)
        {  subdomId = sid;  return;  }

        int get_subdomain_id()
        {  return  subdomId;  }

        void SetDegree(int* deg)
        {
          totnlbf = 1;
          for(int ii=0;ii<DIM;ii++)
          {
            degree[ii] = deg[ii];
            nlbf[ii]   = degree[ii]+1;
            totnlbf   *= nlbf[ii];
          }
          LocalBasisFuncs.resize(totnlbf);
          LocalBasisFuncs.assign(totnlbf, -1);
        }
        
        bool GetFuncFlag()
        {  return FuncFlag; }
        
        void SetFuncFlag(bool flag1)
        {  FuncFlag = flag1; }
        
        void SetKnots(double* knots_)
        {  knots = knots_; }
        
        double* GetKnots(int ind)
        {  return knots[ind]; }

        AABB& GetAABB()
        {  return  bbox; }

        myPoint& GetKnotBegin()
        {  return knotBegin; }

        myPoint& GetKnotEnd()
        {  return knotEnd; }

        myPoint& GetKnotIncrement()
        {  return knotIncr; }

        void SetParent(TreeNode_PTR  parent1)
        {  parent = parent1; }
        
        TreeNode_PTR  GetParent()
        {  return parent; }
        
        int GetNumberOfChildren()
        {  return NUM_CHILDREN; }

        int GetNumberOfNeighbours()
        {  return NUM_NEIGHBOURS; }
        
        void SetNeighbour(int ind, TreeNode_PTR node1)
        {
           assert(ind < pow(2,DIM));
           neighbours[ind] = node1;
        }
        
        TreeNode_PTR  GetNeighbour(int ind) 
        {  return neighbours[ind]; }

        TreeNode_PTR  GetChild(int ind)
        {  return child[ind];	}
        
        static int  GetCount()
        { return nodecount; }

        void SetKnots(int index, double val0, double val1)
        {
          knots[index][0] = val0;
          knots[index][1] = val1;
        }

        void SetKnots(double u0, double u1, double v0, double v1)
        {
          knots[0][0] = u0;
          knots[0][1] = u1;
          knots[1][0] = v0;
          knots[1][1] = v1;
        }

        void SetKnots(double u0, double u1, double v0, double v1, double w0, double w1)
        {
          knots[0][0] = u0;
          knots[0][1] = u1;
          knots[1][0] = v0;
          knots[1][1] = v1;
          knots[2][0] = w0;
          knots[2][1] = w1;
        }

        double  GetKnotSpan(int dir)
        {  return (knots[dir][1] - knots[dir][0]);  }

        double  GetKnotAt(int dir, int loc)
        {  return  knots[dir][loc];  }
        
        bool IsGhost()
        {  return GHOST_FLAG; }
        
        void SetGhostOn()
        {
          GHOST_FLAG  = true;
          ACTIVE_FLAG = false;
        }
        
        void SetGhostOff()
        {  GHOST_FLAG = false;	}
        
        void Activate()
        {  ACTIVE_FLAG = true;	}
        
        void Deactivate()
        {  ACTIVE_FLAG = false; }

        //void SetDomainNumber(int  dd)
        //{ domainNum = dd;  } 

        //int GetDomainNumber()
        //{  return domainNum; }

        int GetDomainNumber()
        {  return ( (domNums.size() == 1) ? domNums[0] : -1); }

        bool IsCutElement()
        {  return (domNums.size() > 1); }

        bool  pointLiesInside(const myPoint& pt);
        
        void  subDivide();
        
        void  unRefine();

        void  plotSelf();
        
        void  printSelf();

        bool  IsBoundary();
        
        bool  IsLeftBoundary();
        bool  IsRightBoundary();
        bool  IsTopBoundary();
        bool  IsBottomBoundary();
        bool  IsFrontBoundary();
        bool  IsBackBoundary();

        void  prepareElemData();

        void  calcSubdivisionMatrix();

        int  computeGaussPointsSubTrias(int ind1=1, int ind2=1, int ind3=1, int ind4=1);

        int  computeGaussPointsAdapIntegration(int ind1=1, int ind2=1, int ind3=1, int ind4=1);

        int  computeGaussPointsAdapIntegrationBoundary(int side, int ind1=1, int ind2=1, int ind3=1, int ind4=1);

        int  checkCutCellValidityAdapIntegration();

        int  prepareCutCell(vector<double>& cutFEMparams);

        int  ResetAdaptiveIntegrationNode();

        int  clearSubtriangulation();

        void  initialiseDOFvalues();
        
        void  checkPartitionOfUnity();

        void  setInitialProfile();

        void  resetMatrixAndVector();

        double  getJacBoundary(int side);

        // for fictitioud domain method
        void  calcStiffnessAndResidualGFEM(MatrixXd& Klocal, VectorXd& Flocal, int domainCur=0);

        void  applyDirichletBCsGFEM(MatrixXd& Klocal, VectorXd& Flocal, int domainCur=0);

        void  applyNeumannBCsGFEM(MatrixXd& Klocal, VectorXd& Flocal, int domainCur=0);

        // for CutFEM for fluid flow 
        void  calcStiffnessAndResidualCutFEMFluid(MatrixXd& Klocal, VectorXd& Flocal, int domainCur=0);

        void  applyDirichletBCsCutFEMFluid(MatrixXd& Klocal, VectorXd& Flocal, int domainCur=0);

        void  applyNeumannBCsCutFEMFluid(MatrixXd& Klocal, VectorXd& Flocal, int domainCur=0);

        void  applyBoundaryConditionsAtApointCutFEMFluid(myDataIntegrateCutFEM&);

        void  applyBoundaryConditionsAtApointCutFEMFluid2(myDataIntegrateCutFEM&);

        //void  setBoundaryGPsCutFEM(int side, vector<double>& boundaryGPs1, vector<double>& boundaryGWs1, double& Jac);
        //void  setBoundaryGPsCutFEM(int side, vector<double>& boundaryGPs1, vector<double>& boundaryGWs1, vector<double>& boundaryGPs2, vector<double>& boundaryGWs2, double& Jac);
        //void  setBoundaryGPsCutFEM(int side, vector<double>& boundaryGPs1, vector<double>& boundaryGWs1, vector<double>& boundaryGPs2, vector<double>& boundaryGWs2, vector<double>& boundaryGPs3, vector<double>& boundaryGWs3, double& Jac);

        // for CutFEM Poisson equation
        void  calcStiffnessAndResidualCutFEMPoisson(MatrixXd& Klocal, VectorXd& Flocal, int domainCur=0);

        void  applyDirichletBCsCutFEMPoisson(MatrixXd& Klocal, VectorXd& Flocal, int domainCur=0);

        void  applyNeumannBCsCutFEMPoisson(MatrixXd& Klocal, VectorXd& Flocal, int domainCur=0);

        void  applyBoundaryConditionsAtApointCutFEMPoisson(myDataIntegrateCutFEM&);

        void  applyBoundaryConditionsAtApointCutFEMPoisson2(myDataIntegrateCutFEM&);

        // for CutFEM ELASTICITY
        void  calcStiffnessAndResidualCutFEMelasticity(MatrixXd& Klocal, VectorXd& Flocal, int domainCur=0);

        void  applyDirichletBCsCutFEMelasticity(MatrixXd& Klocal, VectorXd& Flocal, int domainCur=0);

        void  applyNeumannBCsCutFEMelasticity(MatrixXd& Klocal, VectorXd& Flocal, int domainCur=0);

        void  applyBoundaryConditionsAtApointCutFEMElasticity(myDataIntegrateCutFEM&);

        void  applyBoundaryConditionsAtApointCutFEMElasticity2(myDataIntegrateCutFEM&);

        // for Least-Squares formulation
        void  calcStiffnessAndResidualLSFEM(bool flag, MatrixXd& Klocal, VectorXd& Flocal);

        void  calcResidualLSFEM(VectorXd& Flocal);

        void  applyDirichletBCsLSFEM(bool flag, MatrixXd& Klocal, VectorXd& Flocal);

        void  applyNeumannBCsLSFEM(bool flag, MatrixXd& Klocal, VectorXd& Flocal);

        void  applyBoundaryConditionsAtApoint(myDataIntegrateCutFEM&);

        void  applyGhostPenaltyCutFEM(myDataIntegrateCutFEM&);

        void  computeAndReturnJacobian(int, double*, double*, double*, double, double*, MatrixXd&, VectorXd&, VectorXd&);

        void  AssembleElementVector(int ind, bool flag, double* rhs);
        
        //void  AssembleMatrixAndVector(int, Mat, double*);

        void  AssembleElementMatrix2(int, SparseMatrixXd& globalK);

        void  MatrixToMapResult(int, int, SparseMatrixXd& globalK);

        void  RhsToMapResult(int, int, double*);

        int  calcLoadVector(int ind1=0, int ind2=0, double inp1=0.0, double inp2=0.0);

        int  calcError(int, int domainCur=0);
        
        int  calcErrorPoisson(int, int domainCur=0);

        int  calcErrorElasticity(int, int domainCur=0);

        double  computeTotalBodyForce(int, int);
        double  computeValue(int dir, VectorXd& NN);
        double  computeValuePrev(int dir, VectorXd& NN);
        double  computeValueCur(int dir, VectorXd& NN);
        double  computeValueDot(int dir, VectorXd& NN);
        double  computeValueDotCur(int dir, VectorXd& NN);

        double  computeValueCFM(int domTemp, int dir, VectorXd& NN);

        double  computeValue2(int dir, VectorXd& NN);
        double  computeValue2Prev(int dir, VectorXd& NN);
        double  computeValue2Cur(int dir, VectorXd& NN);

        double  computeForce(int dir, VectorXd& NN);
        double  computeForcePrev(int dir, VectorXd& NN);
        double  computeForceCur(int dir, VectorXd& NN);
        double  computeVorticity(VectorXd& NN);

        void  computeVelocity(myPoint& param, myPoint&  vel);
        void  computeVelocityCur(myPoint& param, myPoint&  vel);

        void  computeVelocityGradient(myPoint& param, MatrixXd&  velGrad);
        void  computeVelocityGradientCur(myPoint& param, MatrixXd&  velGrad);

        void  computeStress(myPoint& param, MatrixXd&  stress);
        void  computeStressCur(myPoint& param, MatrixXd&  stress);

        void  computeVelocityAndStress(myPoint& param, myPoint& vel, MatrixXd&  stress, myPoint& acce);
        void  computeVelocityAndStressCur(myPoint& param, myPoint& vel, MatrixXd&  stress, myPoint& acce);

        void  mapBoundaryPointDataToGlobalBodyForceVector(double* position, double* normal, double arclen, double* useThisData);
};






template <int DIM>
TreeNode<DIM>::TreeNode():FuncFlag(false), level(0), NUM_CHILDREN(0), GHOST_FLAG(false), ACTIVE_FLAG(true), PROCESSED(false)
{
    id = nodecount++;

    child = NULL;
    neighbours = NULL;
    parent = NULL;
    adapIntegNode = NULL;
    
    orientation = ENUM_PARENT;
    
    for(int ii=0;ii<DIM;ii++)
    {
      degree[ii] = 0;
      nlbf[ii] = 0;
    }

    NUM_NEIGHBOURS = 2*DIM;

    neighbours = new TreeNode_PTR[NUM_NEIGHBOURS];
    for(int ii=0;ii<NUM_NEIGHBOURS;ii++)
      neighbours[ii] = NULL;
    
    subdomId = 0;
    //domNums.push_back(0);
}




template <int DIM>
TreeNode<DIM>::TreeNode(int lev):level(lev), NUM_CHILDREN(0), FuncFlag(false), GHOST_FLAG(false), ACTIVE_FLAG(true), PROCESSED(false)
{
    id = nodecount++;

    child = NULL;
    neighbours = NULL;
    parent = NULL;
    adapIntegNode = NULL;

    orientation = ENUM_PARENT;

    for(int ii=0;ii<DIM;ii++)
    {
      degree[ii] = 0;
      nlbf[ii] = 0;
    }

    if(level > MAX_LEVEL)
      MAX_LEVEL = level;

    NUM_NEIGHBOURS = 2*DIM;

    neighbours = new TreeNode_PTR[NUM_NEIGHBOURS];
    for(int ii=0;ii<NUM_NEIGHBOURS;ii++)
      neighbours[ii] = NULL;

    subdomId = 0;
    //domNums.push_back(0);
}


template <int DIM>
TreeNode<DIM>::~TreeNode()
{
  --nodecount;

  if(child != NULL)
  {
    //for(int ii=0;ii<NUM_CHILDREN;ii++)
      //delete child[ii];

    delete [] child;
    child = NULL;
  }

  if(neighbours != NULL)
  {
    //for(int ii=0;ii<NUM_NEIGHBOURS;ii++)
      //delete neighbours[ii];

    delete [] neighbours;
    neighbours = NULL;
  }

  if(parent != NULL)
  {
    //delete parent;
    parent = NULL;
  }

  for(vector<myPoly*>::iterator pObj = subTrias.begin(); pObj != subTrias.end(); ++pObj)
  {
    delete *pObj; // Note that this is deleting what pObj points to, which is a pointer
  }
  subTrias.clear(); // Purge the contents so no one tries to delete them again

}

template <int DIM>
int  TreeNode<DIM>::ResetAdaptiveIntegrationNode()
{
  if(adapIntegNode != NULL)
    delete adapIntegNode;

  adapIntegNode = NULL;
  
  return 1;
}


template <int DIM>
bool  TreeNode<DIM>::IsLeftBoundary()
{
   if(this->IsGhost())
     return false;
   else
     return ((neighbours[LEFT] == NULL) ? false : neighbours[LEFT]->IsGhost() );
}


template <int DIM>
bool  TreeNode<DIM>::IsRightBoundary()
{
   if(this->IsGhost())
     return false;
   else
     return ((neighbours[RIGHT] == NULL) ? false : neighbours[RIGHT]->IsGhost() );
}


template <int DIM>
bool  TreeNode<DIM>::IsTopBoundary()
{
   if(this->IsGhost())
     return false;
   else
     return ((neighbours[NORTH] == NULL) ? false : neighbours[NORTH]->IsGhost() );
}


template <int DIM>
bool  TreeNode<DIM>::IsBottomBoundary()
{
   if(this->IsGhost())
     return false;
   else
     return ((neighbours[SOUTH] == NULL) ? false : neighbours[SOUTH]->IsGhost() );
}


template <int DIM>
bool  TreeNode<DIM>::IsFrontBoundary()
{
   if(this->IsGhost())
     return false;
   else
     return ((neighbours[FRONT] == NULL) ? false : neighbours[FRONT]->IsGhost() );
}


template <int DIM>
bool  TreeNode<DIM>::IsBackBoundary()
{
   if(this->IsGhost())
     return false;
   else
     return ((neighbours[BACK] == NULL) ? false : neighbours[BACK]->IsGhost() );
}


template <int DIM>
void TreeNode<DIM>::resetMatrixAndVector()
{
  //Klocal.setZero();
  //Flocal.setZero();

  return;
}


template<int DIM>
void TreeNode<DIM>::prepareElemData()
{
    int ii;
    totnlbf = 1;
    for(ii=0;ii<DIM;ii++)
    {
      knots[ii][2] = knots[ii][1] - knots[ii][0];
      knots[ii][3] = knots[ii][1] + knots[ii][0];
      
      knotBegin[ii] = knots[ii][0];
      knotEnd[ii]   = knots[ii][1];
      knotIncr[ii]  = knots[ii][2];

      totnlbf *= nlbf[0];
    }

    ndof = GeomData->GetNDOF();
    nsize =  totnlbf * ndof;

    elmDat = &(GeomData->FluidProps[0]);

    //JacMultElem = SolnData->GetJacobianFull() * (0.5*knots[0][2]) * (0.5*knots[1][2]) * (0.5*knots[2][2]);
    //
    // because determinant of the Jacobian is constant for uniform rectangular grids.

    JacMultElem = GeomData->GetJacobianFull();

    for(ii=0;ii<DIM;ii++)
    {
      JacMultElem *= (0.5*knots[ii][2]);

      bbox.minBB[ii] = GeomData->ComputeCoord(ii, knots[ii][0]);
      bbox.maxBB[ii] = GeomData->ComputeCoord(ii, knots[ii][1]);
    }
    //cout << JacMultElem << '\t' << SolnData->GetJacobianFull() << '\t' << (0.5*knots[0][2]) << '\t' << (0.5*knots[1][2]) << '\t' << (0.5*knots[2][2]) << endl;

    return;
}



template<int DIM>
void TreeNode<DIM>::initialiseDOFvalues()
{
    int  ii, jj, ind1, ind2;

    //cout << totnlbf << '\t' << nsize << endl;
    nsize2 = GlobalBasisFuncs.size() * ndof;
    //cout << totnlbf2 << '\t' << nsize2 << endl;

    if(ndof > 1)
    {
      forAssyVec.resize(nsize2);
      for(ii=0;ii<totnlbf2;ii++)
      {
        ind1 = ii * ndof;
        ind2 = GlobalBasisFuncs[ii] * ndof;

        for(jj=0;jj<ndof;jj++)
          forAssyVec[ind1+jj] = ind2 + jj;
      }
      forAssyVec2 = GlobalBasisFuncs;
    }
    else
      forAssyVec = GlobalBasisFuncs;

    counter = 0;

    return;
}



template<int DIM>
void TreeNode<DIM>::AssembleElementMatrix2(int index, SparseMatrixXd& globalMat)
{
  int ii, jj, row, col;

  for(ii=0;ii<nsize;ii++)
  {
    row = forAssyVec[ii];
    for(jj=0;jj<nsize;jj++)
      globalMat.coeffRef(row, forAssyVec[jj]) += Klocal(ii,jj);
  }

  return;
}



template<int DIM>
void TreeNode<DIM>::AssembleElementVector(int ind, bool flag, double* rhs)
{
   for(int ii=0;ii<nsize;ii++)
     rhs[forAssyVec[ii]] += Flocal(ii);

   return;
}



/*
template<int DIM>
void TreeNode<DIM>::AssembleMatrixAndVector(int index, Mat mtx, double* rhs)
{
  PetscErrorCode ierr;
    
  int nn=0, ii, jj;

  //ierr = MatSetValues(mtx,nsize2,&forAssyVec[0],nsize2,&forAssyVec[0],&(Klocal(0,0)),ADD_VALUES);
  // does not work for unsymmetric matrices because of the way Klocal is stored in memory. need to transpose Klocal
  // to use the above code

  for(ii=0;ii<nsize2;ii++)
  {
     rhs[forAssyVec[ii]] += Flocal(ii);
     for(jj=0;jj<nsize2;jj++)
       MatSetValues(mtx, 1, &forAssyVec[ii], 1, &forAssyVec[jj], &(Klocal(ii,jj)), ADD_VALUES);
  }

  return;
}
*/



template<int DIM>
double TreeNode<DIM>::computeValueCFM(int domTemp, int dir, VectorXd& NN)
{
  double  val=0.0;
  for(int ii=0;ii<totnlbf2;ii++)
    val += ( SolnData->solnCFM(domTemp, ndof*GlobalBasisFuncs[ii]+dir) * NN(ii) );

  return val;
}



template<int DIM>
double TreeNode<DIM>::computeValue(int dir, VectorXd& NN)
{
    double  val=0.0;
    double  *Data = &(SolnData->var1(0));

    for(int ii=0;ii<totnlbf2;ii++)
      val += Data[ndof*GlobalBasisFuncs[ii] + dir] * NN(ii);

    return val;
}


template<int DIM>
double TreeNode<DIM>::computeValue2(int dir, VectorXd& NN)
{
    double  val=0.0;
    double  *Data = &(SolnData->var2(0));

    for(int ii=0;ii<totnlbf2;ii++)
      val += Data[GlobalBasisFuncs[ii]] * NN(ii);

    return val;
}



template<int DIM>
double TreeNode<DIM>::computeValueCur(int dir, VectorXd& NN)
{
    double  val=0.0;
    double  *Data = &(SolnData->var1Cur(0));

    for(int ii=0;ii<totnlbf2;ii++)
      val += Data[ndof*GlobalBasisFuncs[ii] + dir] * NN(ii);

    return val;
}



template<int DIM>
double TreeNode<DIM>::computeValue2Cur(int dir, VectorXd& NN)
{
    double  val=0.0;
    double  *Data = &(SolnData->var2Cur(0));

    for(int ii=0;ii<totnlbf2;ii++)
      val += Data[GlobalBasisFuncs[ii]] * NN(ii);

    return val;
}




template<int DIM>
double TreeNode<DIM>::computeValuePrev(int dir, VectorXd& NN)
{
    double  val=0.0;
    double  *Data = &(SolnData->var1Prev(0));

    for(int ii=0;ii<totnlbf2;ii++)
      val += Data[ndof*GlobalBasisFuncs[ii] + dir] * NN(ii);

    return val;
}




template<int DIM>
double TreeNode<DIM>::computeValue2Prev(int dir, VectorXd& NN)
{
    double  val=0.0;
    double  *Data = &(SolnData->var2Prev(0));

    for(int ii=0;ii<totnlbf2;ii++)
      val += Data[GlobalBasisFuncs[ii]] * NN(ii);

    return val;
}




template<int DIM>
double TreeNode<DIM>::computeValueDot(int dir, VectorXd& NN)
{
    double  val=0.0;
    double  *Data = &(SolnData->var1Dot(0));

    for(int ii=0;ii<totnlbf2;ii++)
      val += Data[ndof*GlobalBasisFuncs[ii] + dir] * NN(ii);

    return val;
}


template<int DIM>
double TreeNode<DIM>::computeValueDotCur(int dir, VectorXd& NN)
{
    double  val=0.0;
    double  *Data = &(SolnData->var1DotCur(0));

    for(int ii=0;ii<totnlbf2;ii++)
      val += Data[ndof*GlobalBasisFuncs[ii] + dir] * NN(ii);

    return val;
}


template<int DIM>
double TreeNode<DIM>::computeForce(int dir, VectorXd& NN)
{
    double  val=0.0;
    double  *Data = &(SolnData->bodyForce(0));

    for(int ii=0;ii<totnlbf2;ii++)
      val += Data[ndof*GlobalBasisFuncs[ii] + dir] * NN(ii);

    return val;
}



template<int DIM>
double TreeNode<DIM>::computeForcePrev(int dir, VectorXd& NN)
{
    double  val=0.0;
    double  *Data = &(SolnData->bodyForcePrev(0));

    for(int ii=0;ii<totnlbf2;ii++)
      val += Data[ndof*GlobalBasisFuncs[ii] + dir] * NN(ii);

    return val;
}



template<int DIM>
double TreeNode<DIM>::computeForceCur(int dir, VectorXd& NN)
{
    double  val=0.0;
    double  *Data = &(SolnData->bodyForceCur(0));

    for(int ii=0;ii<totnlbf2;ii++)
      val += Data[ndof*GlobalBasisFuncs[ii] + dir] * NN(ii);

    return val;
}


template<int DIM>
double TreeNode<DIM>::computeVorticity(VectorXd& NN)
{
    double  val=0.0;
    double  *Data = &(SolnData->vorticity(0));

    for(int ii=0;ii<totnlbf2;ii++)
      val += Data[GlobalBasisFuncs[ii]] * NN(ii);

    return val;
}



template<int DIM>
int TreeNode<DIM>::clearSubtriangulation()
{
  for(vector<myPoly*>::iterator pObj = subTrias.begin(); pObj != subTrias.end(); ++pObj)
  {
    delete *pObj; // Note that this is deleting what pObj points to, which is a pointer
  }
  subTrias.clear(); // Purge the contents so no one tries to delete them again

  return 1;
}


#ifdef _DOMAIN1D
typedef TreeNode<1>  node;
#endif


#ifdef _DOMAIN2D
typedef TreeNode<2>  node;
#endif


#ifdef _DOMAIN3D
typedef TreeNode<3>  node;
#endif


#endif
