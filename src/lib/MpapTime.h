
#ifndef incl_MpapTime_h
#define incl_MpapTime_h


#include "MathVector.h"


class MpapTime
{ 
  public:
	
    MpapTime(void) { reset(); return; }

    double cur, prev, prev2, dt, dtPrev, dtMax, write;

    bool dtOK;

    myVector<double> stack;

    void cut(void);
    
    void stck(void);
    
    void reset(void);


// ----- for Deniz's macro163 -------

    myVector<double> Tstep;

    int nTstep;

    void file(char*);

};

#endif

