
#ifndef incl_TimeFunction_h
#define incl_TimeFunction_h


#include "List.h"
#include "TimeFunctionCore.h"


class TimeFunction: public ListItem
{
  public:
	
    TimeFunction(void);
    
    ~TimeFunction();

    int id;
    
    double prop;
    
    List<TimeFunctionCore> fct;

    void update(void);
    
  private:

};


std::ostream &operator<<(std::ostream &, TimeFunction &);

#endif

