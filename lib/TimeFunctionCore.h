
#ifndef incl_TimeFunctionCore_h
#define incl_TimeFunctionCore_h


#include "List.h"


class TimeFunctionCore: public ListItem
{
  public:
	
    TimeFunctionCore(void);
    
    ~TimeFunctionCore();

    double t0, t1, p[10], tp;
    
    double eval(double t);
    
  private:

};

#endif


