
#ifndef incl_SolverMA41_h
#define incl_SolverMA41_h

#include "SolverSparse.h"
#include "SolverWorkSpace.h"


class SolverMA41: public SolverSparse
{	
  public:

    SolverMA41(void);
    
    ~SolverMA41();
    
    double CNTL[10], RINFO[20], *ROWSCA, *COLSCA;

    int    *IS, ICNTL[20], INFO[20], KEEP[50], numSCA, MAXIS;

    virtual int initialise(int p1 = 0, int p2 = 0, int p3 = 0);
   
    virtual int  factorise(void);

    virtual double *solve(double*, int nrhs = 1);

    virtual double *factoriseAndSolve(double*, int nrhs = 1);
    
    virtual void free(void); 
    
  private:

};


#endif



