
#ifndef incl_Element3D10nodedTetrahedron_h
#define incl_Element3D10nodedTetrahedron_h


#include "Element.h"
#include "Element3D.h"
#include "MathVector.h"


class Element3D10nodedTetrahedron: public virtual Element, public Element3D
{
  public:

    Element3D10nodedTetrahedron(void);

    virtual ~Element3D10nodedTetrahedron();
	  
    int ndm(void) { return 3; }
    int nen(void) { return 10; }
    
    virtual void plotOutline(bool defFlg = false);
    
    virtual void paint(bool defFlg = false);

    virtual bool forDomainType(int);

    virtual void putLabel(char*, bool defFlg = false);
    
    virtual void contourPlot(int, int, int, double, double, bool defFlg = true);

    virtual double volume(bool init = false);

    virtual int  nFaces(void) { return 4; }

    virtual void giveFace3D(int, Vector<int> &);

    virtual int  nBasicFacesPerFace(void) { return 4; }

    virtual void defineBasicFace(int, int, int*, unsigned int*);

  private:

	  
};

#endif

