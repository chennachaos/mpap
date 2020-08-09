
#ifndef incl_TimeFunctionCore_h
#define incl_TimeFunctionCore_h

#include <vector>
#include <assert.h>
#include "List.h"


class TimeFunctionCore: public ListItem
{
  public:

    double t0, t1, p[10], tp;

    TimeFunctionCore();

    ~TimeFunctionCore();

    TimeFunctionCore(vector<double>&  inpdata)
    {
      assert(inpdata.size() >= 10);

      t0 = inpdata[0];
      t1 = inpdata[1];

      for(int i=0; i<8; ++i)
        p[i] = inpdata[i+2];
    }

    void  setData(vector<double>&  inpdata);

    double eval(double t);

    double evalValue(double t);

    double evalFirstDerivative(double t);

    double evalSecondDerivative(double t);

  private:

};

#endif


