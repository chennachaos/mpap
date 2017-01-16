
#ifndef incl_MicroCellWulf_h
#define incl_MicroCellWulf_h

#include "Solid.h"


class MicroCellWulf: public Solid
{ 
  public:
    MicroCellWulf(void);                      // constructor

    virtual ~MicroCellWulf();                 // destructor





    
    virtual void prepareInputData(void);
    
    virtual void readInputData(std::ifstream &, MyString &);
    



    virtual void preparePeriodicBoundaryConditions(int);        

    virtual void strainToBoundaryDisplacement(double *);

    virtual void getStressFromReactions(double *);
    
  private:

    int ndPosBnd2D[2];

    int oldNumnp;
};



#include "DomainInlineFunctions.h"

define_reference_cast(microCellWulf,MicroCellWulf)

define_isType(isMicroCellWulf,MICROCELLWULF)
	



	
#endif


