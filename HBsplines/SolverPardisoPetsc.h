
#ifndef incl_SolverPardisoPetsc_h
#define incl_SolverPardisoPetsc_h

#include "SolverPetsc.h"


class SolverPardisoPetsc: public SolverPetsc
{	
  public:

    SolverPardisoPetsc(void);

    ~SolverPardisoPetsc();

    void   *PT[64];

    double DPARM[64], *X, *val;

    PetscInt   *csr, *col;

    PetscScalar *array;

    int    SOLVER, MTYPE, MAXFCT, MNUM, NRHS, MSGLVL, N, IPARM[64], nx;

    virtual int initialise(int p1 = 0, int p2 = 0, int p3 = 0);

    virtual int  factorise();

    virtual double *solve(double*, int nrhs = 1);

    virtual double *factoriseAndSolve(double*, int nrhs = 1);

    virtual int free();

  private:

};


#endif



