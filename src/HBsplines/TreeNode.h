
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

        static int  nodecount, MAX_LEVEL, ndof;

        int  subdomId, level, totnlbf, totnlbf2;
        int  NUM_CHILDREN, NUM_NEIGHBOURS;
        int  id;
        NodeOrientation  orientation;

        double  *elmDat, JacMultElem;

        AABB  bbox;

        myPoint  knotBegin, knotEnd, knotIncr, knotSum;

        bool FuncFlag, GHOST_FLAG, ACTIVE_FLAG, PROCESSED;

        TreeNode_PTR  *child, parent;//, *neighbours;

        TreeNode_PTR  neighbours[2*DIM];

    public:

        vector<int>  domNums, GlobalBasisFuncs, LocalBasisFuncs, LocalBasisFuncsPrev, forAssyVec;
        vector<vector<double> > DirichletData, NeumannData, DerivativeBCData;

        SolutionData  *SolnData;
        GeomDataHBSplines  *GeomData;

        MatrixXd  SubDivMat;

        GaussQuadrature  Quadrature;
        vector<GaussQuadrature>  BoundaryQuadrature;

        vector<myPoly*>  subTrias;

        typedef  AdaptiveOctree<DIM>*  adapTreePTR;

        //typedef  AdaptiveBinarytree<DIM>*  adapTreePTR;

        //AdaptiveOctree<DIM>  *adapIntegNode1;
        //AdaptiveBinarytree<DIM>  *adapIntegNode2;
        adapTreePTR  adapIntegNode;

        TreeNode();

        TreeNode(int lev=0);

        ~TreeNode();

        int  getLevel()
        {  return level; }

        void setLevel(int lev)
        {  level = lev; }

        bool isLeaf()
        {  return (child == NULL); }

        bool isActive()
        {  return  ACTIVE_FLAG; }

        bool isProcessed()
        {  return PROCESSED;}

        int getID()
        {  return id;}

        int getDimension()
        {  return DIM;}

        void  setNdof(int nd)
        { ndof = nd; }

        int  getNdof()
        { return ndof; }

        double  getJacMultElement()
        { return JacMultElem; }

        double  getVolume();

        double  getVolumeGaussPoints(int domTemp);

        void  setSubdomainId(int sid)
        {  subdomId = sid;  return;  }

        int getSubdomainId()
        {  return  subdomId;  }

        int  getLocalBFsSize()
        {  return  totnlbf; }

        void setLocalBFsSize(int dd)
        {
          totnlbf = dd;

          LocalBasisFuncs.resize(totnlbf);
          LocalBasisFuncs.assign(totnlbf, -1);
        }

        void setParent(TreeNode_PTR  parent1)
        {  parent = parent1; }

        TreeNode_PTR  getParent()
        {  return parent; }

        int getNumberOfChildren()
        {  return NUM_CHILDREN; }

        int getNumberOfNeighbours()
        {  return NUM_NEIGHBOURS; }

        void setNeighbour(int ind, TreeNode_PTR node1)
        {
           assert(ind < pow(2,DIM));
           neighbours[ind] = node1;
        }

        TreeNode_PTR  getNeighbour(int ind) 
        {  return neighbours[ind]; }

        TreeNode_PTR  getChild(int ind)
        {  return child[ind];	}

        static int  getCount()
        { return nodecount; }

        AABB& getAABB()
        {  return  bbox; }

        myPoint& getKnotBegin()
        {  return knotBegin; }

        myPoint& getKnotEnd()
        {  return knotEnd; }

        myPoint& getKnotIncrement()
        {  return knotIncr; }

        myPoint& getKnotSum()
        {  return knotSum; }

        void setKnots(int index, double u0, double u1)
        {
          knotBegin[index] = u0;          knotEnd[index]   = u1;

          knotIncr[index] = u1-u0;
          knotSum[index]  = u1+u0;
        }

        void setKnots(double u0, double u1, double v0, double v1)
        {
          knotBegin[0] = u0;          knotEnd[0]   = u1;
          knotIncr[0] = u1-u0;
          knotSum[0]  = u1+u0;

          knotBegin[1] = v0;          knotEnd[1]   = v1;
          knotIncr[1] = v1-v0;
          knotSum[1]  = v1+v0;
        }

        void setKnots(double u0, double u1, double v0, double v1, double w0, double w1)
        {
          knotBegin[0] = u0;          knotEnd[0]   = u1;
          knotIncr[0] = u1-u0;
          knotSum[0]  = u1+u0;

          knotBegin[1] = v0;          knotEnd[1]   = v1;
          knotIncr[1] = v1-v0;
          knotSum[1]  = v1+v0;

          knotBegin[2] = w0;          knotEnd[2]   = w1;
          knotIncr[2] = w1-w0;
          knotSum[2]  = w1+w0;
        }

        double  getKnotspan(int dir)
        {  return  knotIncr[dir];  }

        bool isGhost()
        {  return GHOST_FLAG; }

        void setGhostOn()
        {
          GHOST_FLAG  = true;
          ACTIVE_FLAG = false;
        }

        void setGhostOff()
        {  GHOST_FLAG = false;	}

        void activate()
        {  ACTIVE_FLAG = true;	}

        void deactivate()
        {  ACTIVE_FLAG = false; }

        int getDomainNumber()
        {  return ( (domNums.size() == 1) ? domNums[0] : -1); }

        bool isCutElement()
        {  return (domNums.size() > 1); }

        int getNsize2()
        { return forAssyVec.size(); }

        int getNumberOfBasisFunctions()
        { return GlobalBasisFuncs.size(); }

        int getNumberOfQuadraturePoints()
        {
          if(domNums.size() > 1) // cut-cell
          {
            return  Quadrature.gausspoints.size();
          }
          else
          {
            return  GeomData->gausspoints.size();
          }
        }

        /*
        * Computational effort for each cell.
        * To be used as weight for partitioning the mesh to achieve efficient load-balancing
        */
        int getComputationalEffort()
        {
          if(domNums.size() > 1) // cut-cell
          {
            return  (GlobalBasisFuncs.size() * Quadrature.gausspoints.size());
          }
          else
          {
            return  (GlobalBasisFuncs.size() * GeomData->gausspoints.size());
          }
        }

        bool  pointLiesInside(const myPoint& pt);

        void  subDivide();

        void  unRefine();

        void  plotSelf();

        void  printSelf();

        bool  isBoundary();

        bool  isLeftBoundary();
        bool  isRightBoundary();
        bool  isTopBoundary();
        bool  isBottomBoundary();
        bool  isFrontBoundary();
        bool  isBackBoundary();

        void  prepareElemData();

        void  calcSubdivisionMatrix();

        int  computeGaussPointsSubTrias(int ind1=1, int ind2=1, int ind3=1, int ind4=1);

        int  computeGaussPointsAdapIntegration(int ind1=1, int ind2=1, int ind3=1, int ind4=1);

        int  computeGaussPointsAdapIntegrationBoundary(int side, int ind1=1, int ind2=1, int ind3=1, int ind4=1);

        int  checkCutCellValidityAdapIntegration();

        int  prepareCutCell(vector<double>& cutFEMparams);

        int  resetAdaptiveIntegrationNode();

        int  clearSubtriangulation();

        void  initialiseDOFvalues();

        void  checkPartitionOfUnity();

        void  setInitialProfile();

        double  getJacBoundary(int side);

        // for fictitioud domain method
        void  calcStiffnessAndResidualGFEM(MatrixXd& Klocal, VectorXd& Flocal, int domainCur=0);

        void  applyDirichletBCsGFEM(MatrixXd& Klocal, VectorXd& Flocal, int domainCur=0);

        void  applyNeumannBCsGFEM(MatrixXd& Klocal, VectorXd& Flocal, int domainCur=0);

        void  applyDerivativeBCsGFEM(MatrixXd& Klocal, VectorXd& Flocal, int domainCur=0);

        void  applyBoundaryConditionsAtApoint(myDataIntegrateCutFEM&);

        // for CutFEM for fluid flow 
        void  calcStiffnessAndResidualCutFEMFluid(MatrixXd& Klocal, VectorXd& Flocal, int domainCur=0);

        void  applyDirichletBCsCutFEMFluid(MatrixXd& Klocal, VectorXd& Flocal, int domainCur=0);

        void  applyNeumannBCsCutFEMFluid(MatrixXd& Klocal, VectorXd& Flocal, int domainCur=0);

        void  applyDerivativeBCsCutFEMFluid(MatrixXd& Klocal, VectorXd& Flocal, int domainCur=0);

        void  applyBoundaryConditionsAtApointCutFEMFluid(myDataIntegrateCutFEM&);

        void  applyBoundaryConditionsAtApointCutFEMFluid2(myDataIntegrateCutFEM&);

        void  applyGhostPenaltyCutFEM(myDataIntegrateCutFEM&);

        void  MatrixToMapResult(int, int, SparseMatrixXd& globalK);

        void  RhsToMapResult(int, int, double*);

        int  calcLoadVector(int ind1=0, int ind2=0, double inp1=0.0, double inp2=0.0);

        double  calcError(int, int domainCur=0);

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
TreeNode<DIM>::TreeNode():level(0), NUM_CHILDREN(0), GHOST_FLAG(false), ACTIVE_FLAG(true), PROCESSED(false)
{
    id = nodecount++;

    child = NULL;
    parent = NULL;
    adapIntegNode = NULL;

    orientation = ENUM_PARENT;

    NUM_NEIGHBOURS = 2*DIM;

    //neighbours = new TreeNode_PTR[NUM_NEIGHBOURS];
    for(int ii=0;ii<NUM_NEIGHBOURS;ii++)
      neighbours[ii] = NULL;

    subdomId = 0;
    //domNums.push_back(0);
}




template <int DIM>
TreeNode<DIM>::TreeNode(int lev):level(lev), NUM_CHILDREN(0), GHOST_FLAG(false), ACTIVE_FLAG(true), PROCESSED(false)
{
    id = nodecount++;

    child = NULL;
    parent = NULL;
    adapIntegNode = NULL;

    orientation = ENUM_PARENT;

    if(level > MAX_LEVEL)
      MAX_LEVEL = level;

    NUM_NEIGHBOURS = 2*DIM;

    //neighbours = new TreeNode_PTR[NUM_NEIGHBOURS];
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

  //if(neighbours != NULL)
  //{
    for(int ii=0;ii<NUM_NEIGHBOURS;ii++)
      neighbours[ii] = NULL;
      //delete neighbours[ii];

    //delete [] neighbours;
    //neighbours = NULL;
  //}

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
int  TreeNode<DIM>::resetAdaptiveIntegrationNode()
{
  if(adapIntegNode != NULL)
    delete adapIntegNode;

  adapIntegNode = NULL;

  return 1;
}


template <int DIM>
bool  TreeNode<DIM>::isLeftBoundary()
{
   if(this->isGhost())
     return false;
   else
     return ((neighbours[LEFT] == NULL) ? false : neighbours[LEFT]->isGhost() );
}


template <int DIM>
bool  TreeNode<DIM>::isRightBoundary()
{
   if(this->isGhost())
     return false;
   else
     return ((neighbours[RIGHT] == NULL) ? false : neighbours[RIGHT]->isGhost() );
}


template <int DIM>
bool  TreeNode<DIM>::isTopBoundary()
{
   if(this->isGhost())
     return false;
   else
     return ((neighbours[NORTH] == NULL) ? false : neighbours[NORTH]->isGhost() );
}


template <int DIM>
bool  TreeNode<DIM>::isBottomBoundary()
{
   if(this->isGhost())
     return false;
   else
     return ((neighbours[SOUTH] == NULL) ? false : neighbours[SOUTH]->isGhost() );
}


template <int DIM>
bool  TreeNode<DIM>::isFrontBoundary()
{
   if(this->isGhost())
     return false;
   else
     return ((neighbours[FRONT] == NULL) ? false : neighbours[FRONT]->isGhost() );
}


template <int DIM>
bool  TreeNode<DIM>::isBackBoundary()
{
   if(this->isGhost())
     return false;
   else
     return ((neighbours[BACK] == NULL) ? false : neighbours[BACK]->isGhost() );
}


template<int DIM>
void TreeNode<DIM>::prepareElemData()
{
    elmDat = &(GeomData->FluidProps[0]);

    //JacMultElem = SolnData->getJacobianFull() * (0.5*knotIncr[0]) * (0.5*knotIncr[1]) * (0.5*knotIncr[2]);
    //
    // because determinant of the Jacobian is constant for uniform rectangular grids.

    JacMultElem = GeomData->getJacobianFull();

    for(int ii=0;ii<DIM;ii++)
    {
      JacMultElem *= (0.5*knotIncr[ii]);

      bbox.minBB[ii] = GeomData->computeCoord(ii, knotBegin[ii]);
      bbox.maxBB[ii] = GeomData->computeCoord(ii, knotEnd[ii]);
    }

    return;
}



template<int DIM>
void TreeNode<DIM>::initialiseDOFvalues()
{
    if(ndof > 1)
    {
      int  ii=0, jj=0, ind1=0, ind2=0;
      int  nsize2 = GlobalBasisFuncs.size() * ndof;

      forAssyVec.resize(nsize2);
      for(ii=0;ii<totnlbf2;ii++)
      {
        ind1 = ii * ndof;
        ind2 = GlobalBasisFuncs[ii] * ndof;

        for(jj=0;jj<ndof;jj++)
          forAssyVec[ind1+jj] = ind2 + jj;
      }
    }
    else
      forAssyVec = GlobalBasisFuncs;

    return;
}



/*
template<int DIM>
void TreeNode<DIM>::assembleMatrixAndVector(int index, Mat mtx, double* rhs)
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
  if(adapIntegNode != NULL)
  {
      delete  adapIntegNode;
      adapIntegNode = NULL;
  }

  if( subTrias.size() > 0 )
  {
    for(vector<myPoly*>::iterator pObj = subTrias.begin(); pObj != subTrias.end(); ++pObj)
    {
      delete *pObj; // Note that this is deleting what pObj points to, which is a pointer
    }
    subTrias.clear(); // Purge the contents so no one tries to delete them again
  }

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
