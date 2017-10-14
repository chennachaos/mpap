
#include <cmath>

#include "TimeFunctionCore.h"


TimeFunctionCore::TimeFunctionCore(void)
{
  return;
}



TimeFunctionCore::~TimeFunctionCore()
{
  return;
}



double TimeFunctionCore::eval(double t)
{
  double tol = 1.e-14, tt0, tt1, tt = floor((t+tol-t0)/tp) * tp;

  if (tt < 0.) tt = 0.0;

  tt0 = t0 + tt;
  tt1 = t1 + tt;
  
  if (t < tt0-tol) return 0.0;

  if (t > tt1-tol) return 0.0;

  return ( p[0] + p[1]*(t-tt0) + p[2]*sin(p[3]*(t-tt0)+p[4]) + p[5]*cos(p[6]*(t-tt0)+p[7]) );
}



