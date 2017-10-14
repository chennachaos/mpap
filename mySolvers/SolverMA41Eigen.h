
#ifndef incl_SolverMA41Eigen_h
#define incl_SolverMA41Eigen_h


#include "SolverEigen.h"




class SolverMA41Eigen: public SolverEigen
{
  public:

    SolverMA41Eigen();
    
    ~SolverMA41Eigen();

    vector<int>  row, col;
    vector<double>  array;

    double CNTL[10], RINFO[20], *ROWSCA, *COLSCA;

    int    *IS, ICNTL[20], INFO[20], KEEP[50], numSCA, MAXIS, N, NE;

    virtual int initialise(int p1 = 0, int p2 = 0, int p3 = 0);
   
    virtual int factorise();

    virtual int solve();

    virtual int factoriseAndSolve();

    virtual void free();
    
  private:

};


#endif



