
#ifndef incl_Element3D20nodedBrick_h
#define incl_Element3D20nodedBrick_h


#include "Element.h"
#include "Element3D.h"
#include "MathVector.h"


class Element3D20nodedBrick: public virtual Element, public Element3D
{
  public:

    Element3D20nodedBrick(void);

    virtual ~Element3D20nodedBrick();
	  
    int ndm(void) { return 3; }
    int nen(void) { return 20; }
    
    virtual void plotOutline(bool defFlg = false);
    
    virtual void paint(bool defFlg = false);

    virtual bool forDomainType(int);

    virtual void putLabel(char*, bool defFlg = false);
    
    virtual void contourPlot(int, int, int, double, double, bool defFlg = true);

    virtual double volume(bool init = false);

    virtual int  nFaces(void) { return 6; }

    virtual void giveFace3D(int, Vector<int> &);

    virtual int  nBasicFacesPerFace(void) { return 6; }

    virtual void defineBasicFace(int, int, int*, unsigned int*);

  private:

	  
};

#endif

