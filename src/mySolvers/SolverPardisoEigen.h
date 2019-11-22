
#ifndef incl_SolverPardisoEigen_h
#define incl_SolverPardisoEigen_h

#include "SolverEigen.h"


class SolverPardisoEigen: public SolverEigen
{
  public:

    SolverPardisoEigen(void);

    ~SolverPardisoEigen();

    int   PT[64], IPARM[64];
    int   phase, error, SOLVER, MTYPE, MAXFCT, MNUM, NRHS, MSGLVL;

    double DPARM[64], ddum, *array;

    vector<int>  csr, col, perm;


    virtual int initialise(int p1 = 0, int p2 = 0, int p3 = 0);

    virtual int factorise();

    virtual int solve();

    virtual int factoriseAndSolve();

    virtual int free();
};


#endif



