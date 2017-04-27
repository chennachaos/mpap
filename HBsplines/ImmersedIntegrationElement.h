#ifndef incl_ImmersedIntegrationElement_h
#define incl_ImmersedIntegrationElement_h


#include "util.h"

#include "TreeNode.h"


using namespace std;
using namespace Eigen;

class SolutionData;
class GeomDataLagrange;
class GeomDataHBSplines;



class ImmersedIntegrationElement
{
    private:
    
        static  int  pointcount;

        int DIM, id, nGP;
        
        bool  IS_ACTIVE;

    public:

        vector<int>  pointNums, posIndices, elemNums;
        vector<double>  gausspoints, gaussweights;

        //double  positionOrig[ndim], positionCur[ndim], positionOld[ndim], param[ndim];

        //VectorXd  U, dU, ddU, iU, Un, dUn, ddUn, iUn, force, param;
        myPoint  U, dU, ddU, iU, Un, dUn, ddUn, iUn, force, param;

        //double  normal[ndim], force[ndim], forcePrev[ndim], forceCur[ndim], specVal[ndim];

        VectorXd  Flocal, Flocal2;

        MatrixXd  Khorz, Kvert;

        node*  elem;

        SolutionData  *SolnData;
        GeomDataLagrange  *GeomDataLag;
        GeomDataHBSplines *GeomDataHBS;

        ImmersedIntegrationElement();

        ~ImmersedIntegrationElement();

        int GetID()
        {  return id; }

        static int  GetCount()
        { return pointcount; }

        void SetDimension(int dd)
        {  DIM = dd; return;  }

        int GetDimension()
        {  return DIM;  }

        bool IsActive()
        {  return (IS_ACTIVE == true);  }

        void TurnIsActiveON()
        {  IS_ACTIVE = true;  return ;  }

        void TurnIsActiveOFF()
        {  IS_ACTIVE = false;  return ;  }

        void  printSelf();

        void  findElements();

        void  initialiseDOFvalues();

        void  resetMatrixAndVector();
        
        void  computePointAtGP(int ind, myPoint& pt);

        void  IntegrateForceAndMoment(VectorXd& vectemp, myPoint& centLoc);

        void  IntegrateForceFlexible(int ind1, int ind2, VectorXd& vectemp);

        void  computeKhorzKvertRigid(MatrixXd& tempMat, MatrixXd& tempMat2);

        void  computeKhorzKvertFlexible(int ind1, int ind2, MatrixXd& tempMat, MatrixXd& tempMat2);

        void  calcStiffnessAndResidual(int ind1=0, int ind2=0, double inp1=0.0, double inp2=0.0);

        void  AssembleElementMatrix(int, SparseMatrixXd&);

        void  AssembleElementVector(int ind, bool flag, double* rhs);

        void  AssembleMatrixAndVector(int, SparseMatrixXd&, double*);

        void  prepareElemData();

        void  computeBodyForce(bool, double*);

        void  getDataFromGlobalSolutionVector(bool flag);

        void  mapDataToGlobalBodyForceVector(bool flag, double* useThisData);

        void reset();

        void  computeVelocity(const VectorXd& NN, myPoint&  velSpec);

        void  computeVelocityCur(const VectorXd& NN, myPoint&  velSpec);

        void  computeAcceleration(const VectorXd& NN, myPoint&  velSpec);

        void  computeAccelerationCur(const VectorXd& NN, myPoint&  velSpec);
};






#endif


