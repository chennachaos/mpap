
#ifndef incl_SolverPardisoPetsc_h
#define incl_SolverPardisoPetsc_h

#include "SolverPetsc.h"


class SolverPardisoPetsc: public SolverPetsc
{	
  public:

    SolverPardisoPetsc(void);

    ~SolverPardisoPetsc();

    void   *PT[64];

    double DPARM[64], *val, ddum;
    vector<double>  solnX;

    PetscInt   *csr, *col, *perm;

    PetscScalar *array;

    int    phase, error, SOLVER, MTYPE, MAXFCT, MNUM, NRHS, MSGLVL, IPARM[64];

    virtual int initialise(int p1 = 0, int p2 = 0, int p3 = 0);

    virtual int  factorise();

    virtual int  solve();

    virtual int  factoriseAndSolve();

    virtual int free();

  private:

};


#endif



