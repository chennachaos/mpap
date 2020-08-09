
#ifndef incl_SolverPardisoPetsc_h
#define incl_SolverPardisoPetsc_h

#include "SolverPetsc.h"


class SolverPardisoPetsc: public SolverPetsc
{	
  public:

    SolverPardisoPetsc(void);

    ~SolverPardisoPetsc();

    int   PT[64], IPARM[64];
    int   phase, error, SOLVER, MTYPE, MAXFCT, MNUM, NRHS, MSGLVL;

    double DPARM[64], ddum;

    vector<double>  rhsTemp, solnTemp;

    PetscInt   *csr, *col, *perm;

    PetscScalar *array;

    virtual int initialise(int p1 = 0, int p2 = 0, int p3 = 0);

    virtual int  factorise();

    virtual int  solve();

    virtual int  factoriseAndSolve();

    virtual int free();

  private:

};


#endif



