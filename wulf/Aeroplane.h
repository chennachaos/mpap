
#ifndef incl_Aeroplane_h
#define incl_Aeroplane_h


#include "Domain.h"
#include "DAWindow.h"
#include "DAFunctionPlot.h"


class Aeroplane: public Domain
{ 
  public:

    Aeroplane(void);
    virtual ~Aeroplane();

    DAWindow *aDAW;

    DAFunctionPlot *flightProfile;

    virtual void readInputData(std::ifstream &, MyString &);

    virtual void prepareInputData(void);

    virtual void prepareInteractions(void);

    virtual void initFlightSimulatorDisplay(bool*);

    virtual int factoriseSolveAndUpdate(void);

  private:


};




#include "DomainInlineFunctions.h"

define_reference_cast(aeroplane,Aeroplane)

define_isType(isAeroplane,AEROPLANE)
	



#endif




