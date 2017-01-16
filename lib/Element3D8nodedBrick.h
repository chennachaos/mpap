
#ifndef incl_Element3D8nodedBrick_h
#define incl_Element3D8nodedBrick_h


#include "Element.h"
#include "Element3D.h"
#include "MathVector.h"


class Element3D8nodedBrick: public virtual Element, public Element3D
{
  public:

    Element3D8nodedBrick(void);

    virtual ~Element3D8nodedBrick();
	  
    int ndm(void) { return 3; }
    int nen(void) { return 8; }
    
    virtual void plotOutline(bool defFlg = false);
    
    virtual void paint(bool defFlg = false);

    virtual bool forDomainType(int);

    virtual void putLabel(char*, bool defFlg = false);
    
    virtual void contourPlot(int, int, int, double, double, bool defFlg = true);

    virtual double volume(bool init = false);

    virtual int  nFaces(void) { return 6; }

    virtual void giveFace3D(int, Vector<int> &);

    virtual int  nBasicFacesPerFace(void) { return 2; }

    virtual void defineBasicFace(int, int, int*, unsigned int*);

  private:

	  
};

#endif

