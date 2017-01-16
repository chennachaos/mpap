
#ifndef incl_VortexSheet_h
#define incl_VortexSheet_h


#include "Domain.h"
#include "MathMatrix.h"

using namespace  MatricesWulf;

class VortexSheet: public Domain
{ 
  public:

    VortexSheet(void);
    virtual ~VortexSheet();

    virtual void readInputData(std::ifstream &, MyString &);

    virtual void prepareInputData(void);

    virtual void prepareInteractions(void);

    virtual void doForVortexSheet(void);

  private:

    MatrixFullArray<double> x;

    double v, rho, alpha;

    int nSheet, iSep;

};




#include "DomainInlineFunctions.h"

define_reference_cast(vortexSheet,VortexSheet)

define_isType(isVortexSheet,VORTEXSHEET)
	



#endif




