
#ifndef incl_MicroCell_h
#define incl_MicroCell_h

#include "Solid.h"

enum {linearbc, tractionbc, periodicbc, clinearbc, ctractionbc};


class MicroCell: public Solid
{ 
  public:
    MicroCell(void);                      // constructor

	MicroCell(const MicroCell& p);         // copy constructor

        virtual ~MicroCell();                 // destructor

	virtual void prepareInputData(void);

    virtual void prepareInteractions(void);
    
    virtual void readInputData(std::ifstream &, MyString &);





  private:


};



#include "DomainInlineFunctions.h"

define_reference_cast(microCell,MicroCell)

define_isType(isMicroCell,MICROCELL)
	



	
#endif


