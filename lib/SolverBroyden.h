
#ifndef incl_SolverBroyden_h
#define incl_SolverBroyden_h

#include "Solver.h"
#include "SolverWorkSpace.h"
#include "List.h"
#include "MathVector.h"
#include "MathMatrix.h"
#include "DependentDoF.h"
#include "Element.h"


class SolverBroyden: public Solver
{	
  public:

    SolverBroyden(void);
    
    ~SolverBroyden();

    virtual int  initialise(int p1 = 0, int p2 = 0, int p3 = 0);
   
    virtual double *factoriseAndSolve(double*, int);

    virtual void free(void);
    
  private:

    int neq;

    double *mtx, *dxn, *dFn, *Fn;
};


#endif



