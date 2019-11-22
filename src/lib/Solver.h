
#ifndef incl_Solver_h
#define incl_Solver_h

#include <iostream>

#include "MathVector.h"
#include "MathMatrix.h"
#include "List.h"
#include "Domain.h"


using namespace std;




class Solver
{
  public:

    Solver(void);
	  
    virtual ~Solver();

    int  currentStatus, whatToSolveFor;

    bool checkIO;

    Domain *dom;  // pointer to allow fast access to associated domain data
    
    virtual int initialise(int p1 = 0, int p2 = 0, int p3 = 0) = 0;

    virtual bool isChildOfSolverSparse(void) { return false; }

    virtual void prepareMatrixPattern(Domain *)
      { cout << "  'prepareMatrixPattern' is not available for this solver!\n\n"; return; }

    virtual void zeroMtx(void)
      { cout << "  'zeroMtx' is not available for this solver!\n\n"; return; }

    //virtual int assembleElemMtx(Element*)
      //{ cout << "  'assembleElemMtx' is not available for this solver!\n\n"; return -1; } 
    
    //virtual int assembleElemVec(Element*, bool)
      //{ cout << "  'assembleElemVec' is not available for this solver!\n\n"; return -1; } 
    
    virtual int  factorise(void)
      { cout << "  'factorise' is not available for this solver!\n\n"; return -1; } 

    virtual double *solve(double*, int nrhs = 1)
      { cout << "  'solve' is not available for this solver!\n\n"; return NULL; } 

    virtual double *factoriseAndSolve(double*, int nrhs = 1)
      { cout << "  'factoriseAndSolve' is not available for this solver!\n\n"; return NULL; } 

    virtual void free(void)
      { cout << "  'free' is not available for this solver!\n\n"; return; } 

    virtual void printInfo(void)
      { cout << "  'printInfo' is not available for this solver!\n\n"; return; }

    virtual void printMatrix(int,int,bool,int indent = 0, bool interactive = false)
      { cout << "  'printMatrix' is not available for this solver!\n\n"; return; }
    
    virtual double giveMatrixCoefficient(int,int)
      { cout << "  'giveMatrixCoefficient' is not available for this solver!\n\n"; return 0.; }

    virtual void copyToSimpleMatrix(double*&)
      { cout << "  'copyToSimpleMatrix' is not available for this solver!\n\n"; return; }

  private:

};



enum { EMPTY, PATTERN_OK, INIT_OK, ASSEMBLY_OK, FACTORISE_OK };


#endif



