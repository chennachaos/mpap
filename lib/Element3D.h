
#ifndef incl_Element3D_h
#define incl_Element3D_h


#include "Element.h"
#include "MathVector.h"



class Element3D: public virtual Element
{
  public:

    //Element3D(void);

    //virtual ~Element3D();

    virtual int nFaces(void) = 0;

    virtual void giveFace3D(int, Vector<int> &) = 0;

    virtual int  nBasicFacesPerFace(void) = 0;

    virtual void defineBasicFace(int, int, int*, unsigned int*) = 0;

  private:

	  
};

#endif

