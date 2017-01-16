
#ifndef incl_FlexibleWing_h
#define incl_FlexibleWing_h


#include "Domain.h"
#include "MathMatrix.h"
#include "MyString.h"
#include "SolverMA41.h"


class FlexibleWing: public Domain
{ 
  public:

    FlexibleWing(void);

    virtual ~FlexibleWing();

    virtual void readInputData(std::ifstream &, MyString &);

    virtual void prepareInputData(void);

    virtual void prepareInteractions(void);

    virtual void doForFlexibleWing(int, bool, double *, MyString &);

  private:

    VectorArray<double> y, c, da, e, GJ, l, 
                        A, B, f1, f2, w, r0, sol, ra0, dsol, dCL, dCDi, dCD, ddCDi, ccl, ccdi, ccd;

    MatrixSparse<double> rda;

    VectorArray<int> aerofoil;

    int numsc;

    double S, q, *ai, *th;

    SolverMA41 *solver;

    void prepareForLiftingLine(void);

    void prepareForFSI(void);

    void generateSystemMatrixAndPrepareRHS(bool);

    void generateFullRHS(double *, double);

    void solveForAlpha0(double &, double);

    void calculateLiftAndDrag(double &, double &, double &);

    void printResults(double, double, double, double);

    void writeResults(MyString &);

    void calculateDerivatives0(void);

    void calculateDerivatives1(void);

    void diffTest(double);

    void optimiseTwist(bool, double, double &, double);
};




#include "DomainInlineFunctions.h"

define_reference_cast(flexibleWing,FlexibleWing)

define_isType(isFlexibleWing,FLEXIBLEWING)
	



#endif




