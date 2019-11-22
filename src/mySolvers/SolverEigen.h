
#ifndef incl_SolverEigen_h
#define incl_SolverEigen_h


#include "Solver.h"
#include "headersEigen.h"
#include <vector>
#include "myIncompleteLUT.h"

using namespace myIterSolvers;

enum { PARDISO_STRUCT_SYM, PARDISO_UNSYM };


class SolverEigen
{
  public:

    VectorXd   rhsVec, soln, solnPrev, var1, var1Prev, var2, var2Prev, var3, var3Prev;
    VectorXd   f1, f2, f3, reac;

    SparseMatrixXd  mtx, matA, matB, matC, matE, PreCondSchur;

    int nRow, nCol, nnz, currentStatus, algoType, nU, nP, nL, nS, update_precond;

    bool checkIO, STABILISED;

    double normPrev, normCur, normRef; // norm of solution error

    myIncompleteLUT<double>  precond;

    ////////////////////////////
    //
    // member functions
    //
    ///////////////////////////

    SolverEigen();

    virtual ~SolverEigen();

    virtual int initialise(int p1 = 0, int p2 = 0, int p3 = 0);

    int setSolverAndParameters();

    void setAlgorithmType(int tt)
    {  algoType = tt; return; }

    virtual void printMatrixPatternToFile();

    virtual void zeroMtx();

    virtual int free();

    virtual void printInfo();

    virtual void printMatrix(int dig=8, int dig2=4, bool gfrmt=true, int indent = 0, bool interactive = false);

    virtual double giveMatrixCoefficient(int,int);

    virtual int assembleMatrixAndVector(vector<int>& row, vector<int>& col, MatrixXd& Klocal, VectorXd& Flocal);

    virtual int assembleMatrixAndVector(int r1, int c1, vector<int>& row, vector<int>& col, MatrixXd& Klocal, VectorXd& Flocal);

    virtual int assembleMatrixAndVector(int r1, int c1, vector<int>& row, MatrixXd& Klocal, VectorXd& Flocal);

    virtual int assembleMatrixAndVectorCutFEM(int r1, int c1, vector<int>& tempVec, vector<int>& forAssy, MatrixXd& Klocal, VectorXd& Flocal);

    virtual int assembleMatrixAndVectorCutFEM2(int r1, int c1, vector<int>& tempVec, 
                vector<int>& forAssy1, vector<int>& forAssy2, MatrixXd& Klocal, VectorXd& Flocal1, VectorXd& Flocal2);

    virtual int assembleMatrixAndVectorCutFEM3(int r1, int c1, vector<int>& row, vector<int>& col,
                vector<int>& forAssy1, vector<int>& forAssy2, MatrixXd& Klocal, VectorXd& Flocal1);

    virtual int assembleVector(int r1, int c1, vector<int>& row, VectorXd& Flocal);

    virtual int factorise();

    virtual int solve();

    virtual int factoriseAndSolve();

    void setupMatricesAndVectors();

    void solverSchurCG();

    void solverSchurGMRES();

    void solverSchurBiCGSTAB();

    void solverUzawaType1();

    void solverUzawaType2();

    int  myBiCGSTAB();

    void resetPrecondFlag()
    {
      update_precond = 1;
    }

    void  updatePreconditioner();

    void  computeConditionNumber();

};



#endif

