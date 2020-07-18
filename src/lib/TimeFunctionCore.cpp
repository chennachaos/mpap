
#include <cmath>
#include <assert.h>

#include "TimeFunctionCore.h"


TimeFunctionCore::TimeFunctionCore()
{
  return;
}


TimeFunctionCore::~TimeFunctionCore()
{
  return;
}


void  TimeFunctionCore::setData(std::vector<double>&  inpdata)
{
  assert(inpdata.size() >= 10);

  t0 = inpdata[0];
  t1 = inpdata[1];

  for(int i=0; i<8; ++i)
    p[i] = inpdata[i+2];

  return;
}



double TimeFunctionCore::eval(double t)
{
  double tol = 1.e-14;

  if (t < t0-tol) return 0.0;
  if (t > t1-tol) return 0.0;

  return ( p[0] + p[1]*(t-t0) + p[2]*sin(p[3]*(t-t0)+p[4]) + p[5]*cos(p[6]*(t-t0)+p[7]) );
}



double TimeFunctionCore::evalValue(double t)
{
  return eval(t);
}



double TimeFunctionCore::evalFirstDerivative(double t)
{
  double tol = 1.e-14;

  if (t < t0-tol) return 0.0;
  if (t > t1-tol) return 0.0;

  return ( p[1] + p[2]*p[3]*cos(p[3]*(t-t0)+p[4]) - p[5]*p[6]*sin(p[6]*(t-t0)+p[7]) );
}



double TimeFunctionCore::evalSecondDerivative(double t)
{
  double tol = 1.e-14;

  if (t < t0-tol) return 0.0;
  if (t > t1-tol) return 0.0;

  return ( -p[2]*p[3]*p[3]*sin(p[3]*(t-t0)+p[4]) - p[5]*p[6]*p[6]*cos(p[6]*(t-t0)+p[7]) );
}




