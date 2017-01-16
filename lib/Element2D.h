
#ifndef incl_Element2D_h
#define incl_Element2D_h


#include "Element.h"
#include "MathVector.h"



class Element2D: public virtual Element
{
  public:

    //Element2D(void);

    //virtual ~Element2D();

    virtual void givePlotSequence2D(Vector<int> &) = 0;

  private:

	  
};

#endif

