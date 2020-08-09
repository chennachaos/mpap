
#include <iostream>

#include "TimeFunction.h"
#include "MpapTime.h"


extern MpapTime mpapTime;


TimeFunction::TimeFunction()
{
  id = 0;

  return;
}



TimeFunction::~TimeFunction()
{
  return;
}



void TimeFunction::update()
{
  prop = 0.;
  for(int i=0; i<fct.n; i++)  prop += fct[i].eval(mpapTime.cur);

  return;
}



// overload printing to screen

std::ostream &operator<<(std::ostream &os, TimeFunction &tmFct)
{
  if (tmFct.fct.n == 0) return os;

  double *p;

  for(int i=0; i<tmFct.fct.n; i++)
  {
    p = tmFct.fct[i].p;

    os << " id = " << tmFct.id << 
      ": " << tmFct.fct[i].t0 << " <= t < " << tmFct.fct[i].t1 <<
      " :  lam[" << i+1 << "](t) = " << p[0] << " + " << p[1] << " * t + " 
             << p[2] << " * sin(" << p[3] << "*t+" << p[4] << " + "
             << p[5] << " * cos(" << p[6] << "*t+" << p[7] << ")\n";
  }

  return os;
}





