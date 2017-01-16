
#ifndef incl_DAWindow_h
#define incl_DAWindow_h


#include "List.h"


class DAWindow: public ListItem
{
  public:

    DAWindow(void);

    ~DAWindow();
  
    void setup(int, int, int, int, double, double, bool keepAspectRatio = false);

    void wipe(int colour = 8);

    void frame(int colour = 7, int d = 1);

    void fillRectangle(int, double, double, double, double, bool flag = false);

    int x0, y0, dx, dy;

    double mx, my;

  private:


};




#endif

