#ifndef incl_ContactElementPointToPoint2D_h
#define incl_ContactElementPointToPoint2D_h

#include "SolverEigen.h"

#include <vector>

using std::vector;

class  ImmersedSolid;
class  SolutionData;
//class  SolutionDataFluid ;
//class  SolutionDataSolid ;



class ContactElementPointToPoint2D
{
    private:
    
        //static  int  pointcount;

        int DIM, id, ndof, nNode, totalDOF, size;

    public:

        vector<int>  nodeNumsLocal, nodeNumsGlobal, posIndices, assy4r;
        vector<vector<int> >  boundaryConds;
        
        vector<ImmersedSolid*> rigidbodypointers;

        vector<vector<double> >  nodePos;

        VectorXd  Flocal;

        MatrixXd  Klocal;

        SolutionData  *FluidSolnData;
        SolutionData  *SolidSolnData;

        ContactElementPointToPoint2D();

        ~ContactElementPointToPoint2D();

        int getID()
        {  return id; }

        //static int  getCount()
        //{ return pointcount; }

        int getDimension()
        {  return DIM;  }

        void  printSelf();

        void  findElements();

        void  initialiseDOFvalues();

        void  calcStiffnessAndResidual(int ind1=0, int ind2=0, double inp1=0.0, double inp2=0.0);

        void  assembleElementMatrix(int, SparseMatrixXd&);

        void  assembleElementVector(int ind, bool flag, double* rhs);

        void  AssembleMatrixAndVector(int, SparseMatrixXd&, double*);

        void  prepareElemData();

        void reset();
};






#endif


