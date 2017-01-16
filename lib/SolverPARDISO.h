
#ifndef incl_SolverPARDISO_h
#define incl_SolverPARDISO_h

#include "SolverSparse.h"
#include "SolverWorkSpace.h"


class SolverPARDISO: public SolverSparse
{	
  public:

    SolverPARDISO(void);
    
    ~SolverPARDISO();

    void   *PT[64];

    double DPARM[64], *X;

    int    SOLVER, MTYPE, MAXFCT, MNUM, NRHS, MSGLVL, N, IPARM[64], nx;

    virtual int initialise(int p1 = 0, int p2 = 0, int p3 = 0);
   
    virtual int  factorise(void);

    virtual double *solve(double*, int nrhs = 1);

    virtual double *factoriseAndSolve(double*, int nrhs = 1);
    
    virtual void free(void); 
    
  private:

};



enum { PARDISO_STRUCT_SYM, PARDISO_UNSYM };


#endif



