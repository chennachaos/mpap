
#ifndef incl_LiftingLine_h
#define incl_LiftingLine_h


#include "Domain.h"
#include "MathMatrix.h"


class LiftingLine: public Domain
{ 
  public:

    LiftingLine(void);

    virtual ~LiftingLine();

    virtual void readInputData(std::ifstream &, MyString &);

    virtual void prepareInputData(void);

    virtual void prepareInteractions(void);

    virtual void doForLiftingLine(void);

  private:

    VectorArray<double> y, c, da, l;

    VectorArray<int> aerofoil;

    int numsc;

    double a0, S;

    void analyseWing(double *, double *, double *, int *,
                     double *, double *, double *,
                     double &, double &, double &);

    void calculateMatrixAndRHS(double *, double *);

    void calculateLiftAndDrag(double *, double *, double *, double *,
                              double &, double &, double &);

    void printResults(double *, double *, double *, double &, double &, double &);

    void plotResults(double *, double *, double *, double &, double &, double &);

    void calculateDerivatives(double *, double *, double *, int *, double *, double *,
                              double *, double *, double *, double *, double *);

    void diffTestHlp(double *, double *, double *, int *, double *, double *, double *, 
                     double *, double *);

    void diffTest(double *, double *, double *, int *, double *, double *, double *);

    void optimiseTwist(double *, double *, double *, int *, double *, double *, double *);
};




#include "DomainInlineFunctions.h"

define_reference_cast(liftingLine,LiftingLine)

define_isType(isLiftingLine,LIFTINGLINE)
	



#endif




