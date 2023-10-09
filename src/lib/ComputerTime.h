
#ifndef incl_ComputerTime_h
#define incl_ComputerTime_h


#include "MyStringList.h"
#include "MathVector.h"



class ComputerTime
{
  public:

    ComputerTime(void);
    
    ~ComputerTime();

    void free(void);

    void reset(void);
    
    void sleep(double ms);

    bool go(char*, bool warn = true);

    double stop(char*);

    double stopAndPrint(char*);

  private:

    MyStringList name;

    myVector<unsigned int> t0;

    double clocksPerSec, secPerClock;
};

    


#endif

