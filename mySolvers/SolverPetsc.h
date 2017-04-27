
#ifndef incl_SolverPetsc_h
#define incl_SolverPetsc_h

//#include "headersBasic.h"

#include "SolverEigen.h"

#include "Solver.h"

#include "petscksp.h"
#include "petscmat.h"

//#include <Eigen/Dense>

using namespace std;
//using Eigen::VectorXd;
//using Eigen::MatrixXd;

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
    
    int SetSolverAndParameters();

    virtual bool isChildOfSolverSparse(void) { return false; }

    virtual int zeroMtx();

    virtual int free();

    virtual int printInfo();

    virtual int printMatrix(int dig=8, int dig2=4, bool gfrmt=true, int indent = 0, bool interactive = false);

    virtual double giveMatrixCoefficient(int,int);

    virtual int AssembleMatrixAndVector(vector<int>& row, vector<int>& col, MatrixXd& Klocal, VectorXd& Flocal);

    virtual int AssembleMatrixAndVector(int r1, int c1, vector<int>& row, vector<int>& col, MatrixXd& Klocal, VectorXd& Flocal);

    virtual int AssembleMatrixAndVectorCutFEM(int r1, int c1, vector<int>& tempVec, vector<int>& forAssy, MatrixXd& Klocal, VectorXd& Flocal);

    virtual int AssembleMatrixAndVectorCutFEM2(int r1, int c1, vector<int>& tempVec, 
                vector<int>& forAssy1, vector<int>& forAssy2, MatrixXd& Klocal, VectorXd& Flocal1, VectorXd& Flocal2);

    virtual int AssembleMatrixAndVectorCutFEM3(int r1, int c1, vector<int>& row, vector<int>& col,
                vector<int>& forAssy1, vector<int>& forAssy2, MatrixXd& Klocal, VectorXd& Flocal1);

    virtual int AssembleVector(int r1, int c1, vector<int>& row, VectorXd& Flocal);

    virtual int AssembleMatrixAndVectorMixedFormulation(int r1, int c1, vector<int>& vec1, vector<int>& vec2, MatrixXd& Klocal, VectorXd& Flocal);

    //virtual double *solve(double*, int nrhs = 1);
    //virtual double *factoriseAndSolve(double*, int nrhs = 1);

    virtual int  factorise();

    virtual int solve();

    virtual int factoriseAndSolve();
    
    int  solveSerial(SparseMatrixXd& matEigen, VectorXd& rhsPetsc, VectorXd& solnPetsc);
    
    int  solveParallel(SparseMatrixXd& matEigen, VectorXd& rhsPetsc, VectorXd& solnPetsc);

};



#endif

