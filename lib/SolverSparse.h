
#ifndef incl_SolverSparse_h
#define incl_SolverSparse_h

#include "Solver.h"
#include "SolverWorkSpace.h"
#include "List.h"
#include "MathVector.h"
#include "MathMatrix.h"
#include "DependentDoF.h"
#include "Element.h"

using namespace MatricesWulf;

class SolverSparse: public Solver
{	
  public:

    SolverSparse(void);
    
    ~SolverSparse();

    MatrixSparseArray<double> mtx;

    bool comprMtxFlg;

    VectorArray<int> compr;

    virtual bool isChildOfSolverSparse(void) { return true; }

    virtual void prepareMatrixPattern(Domain *);

    virtual void zeroMtx(void);
    
    virtual int assembleElemMtx(Element*);
    
    virtual int assembleElemVec(Element*, bool);

    virtual void free(void);

    virtual void printInfo(void);

    virtual void printMatrix(int,int,bool,int indent = 0, bool interactive = false);

    virtual double giveMatrixCoefficient(int,int);
    
    virtual void copyToSimpleMatrix(double*&);

  private:

    bool uDepPosExists(Element **,
                       VectorArray<int> *,
                       ListArray< List< Vector<int> > > &,
                       List< Vector<int> > *, 
		       List< Vector<int> > *,
		       DependentDoF *, 
		       int, int, int, bool, int*);
    
};


#endif



