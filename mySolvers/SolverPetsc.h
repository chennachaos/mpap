
#ifndef incl_SolverPetsc_h
#define incl_SolverPetsc_h

#include "SolverEigen.h"
#include "Solver.h"
#include "petscksp.h"
#include "petscmat.h"


using namespace std;


class SolverPetsc
{
  public:

    Vec  rhsVec, soln, solnPrev, reac;
    Mat  mtx; // linear system matrix
    KSP  ksp; // linear solver context
    PC   pc; // preconditioner context

    PetscInt nRow, nCol, nnz;
    
    MatInfo info;

    int currentStatus;

    bool  checkIO;

    PetscReal norm; // norm of solution error
    
    PetscErrorCode ierr;

    PetscMPIInt size;
    
    PetscViewer    viewer_matx, viewer_vect;

    ////////////////////////////
    //
    // member functions
    //
    ///////////////////////////
    

    SolverPetsc();

    virtual ~SolverPetsc();

    virtual int initialise(int p1 = 0, int p2 = 0, int p3 = 0);
    
    int setSolverAndParameters();

    virtual bool isChildOfSolverSparse(void) { return false; }

    virtual int zeroMtx();

    virtual int free();

    virtual int printInfo();

    virtual int printMatrix(int dig=8, int dig2=4, bool gfrmt=true, int indent = 0, bool interactive = false);

    virtual double giveMatrixCoefficient(int,int);

    virtual int assembleMatrixAndVector(vector<int>& row, vector<int>& col, MatrixXd& Klocal, VectorXd& Flocal);

    virtual int assembleMatrixAndVector(int r1, int c1, vector<int>& row, vector<int>& col, MatrixXd& Klocal, VectorXd& Flocal);

    virtual int assembleMatrixAndVectorCutFEM(int r1, int c1, vector<int>& tempVec, vector<int>& forAssy, MatrixXd& Klocal, VectorXd& Flocal);

    virtual int assembleMatrixAndVectorCutFEM2(int r1, int c1, vector<int>& tempVec, 
                vector<int>& forAssy1, vector<int>& forAssy2, MatrixXd& Klocal, VectorXd& Flocal1, VectorXd& Flocal2);

    virtual int assembleMatrixAndVectorCutFEM3(int r1, int c1, vector<int>& row, vector<int>& col,
                vector<int>& forAssy1, vector<int>& forAssy2, MatrixXd& Klocal, VectorXd& Flocal1);

    virtual int assembleVector(int r1, int c1, vector<int>& row, VectorXd& Flocal);

    virtual int assembleMatrixAndVectorMixedFormulation(int r1, int c1, vector<int>& vec1, vector<int>& vec2, MatrixXd& Klocal, VectorXd& Flocal);

    virtual int factorise();

    virtual int solve();

    virtual int factoriseAndSolve();
    
    int  solveSerial(SparseMatrixXd& matEigen, VectorXd& rhsPetsc, VectorXd& solnPetsc);
    
    int  solveParallel(SparseMatrixXd& matEigen, VectorXd& rhsPetsc, VectorXd& solnPetsc);

};



#endif

