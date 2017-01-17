
#ifndef incl_SolverPardisoEigen_h
#define incl_SolverPardisoEigen_h

#include "SolverEigen.h"


class SolverPardisoEigen: public SolverEigen
{
  public:

    SolverPardisoEigen(void);

    ~SolverPardisoEigen();

    void   *PT[64];

    double  DPARM[64], *array;

    vector<int>  csr, col;
    vector<int> perm;

    int    SOLVER, MTYPE, MAXFCT, MNUM, NRHS, MSGLVL, IPARM[64];

    virtual int initialise(int p1 = 0, int p2 = 0, int p3 = 0);

    virtual int factorise();

    virtual int solve();

    virtual int factoriseAndSolve();

    virtual void free();

};


#endif



