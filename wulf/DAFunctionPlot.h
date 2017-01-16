
#ifndef incl_DAFunctionPlot_h
#define incl_DAFunctionPlot_h


#include "DAWindow.h"
#include "MathVector.h"


class DAFunctionPlot: public DAWindow
{
  public:

    DAFunctionPlot(void);

    ~DAFunctionPlot();

    virtual void setup(int, int, int, int, 
                       int fc = 0, int bc = 8, int gc = 7, int maxgl = 10, 
                       bool kARatio = false, bool grid = true);

    void draw(void);

    void addXY(double, double);    

    Vector<double> x, y;

  private:

    double mnx, mny,
           mxx, mxy;

    int mxn, bCol, gCol, fCol, mxgl;

    bool kAR, showGrid;  // this remains to be implemented!!

};




#endif

