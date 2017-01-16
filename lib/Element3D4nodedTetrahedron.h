
#ifndef incl_Element3D4nodedTetrahedron_h
#define incl_Element3D4nodedTetrahedron_h


#include "Element.h"
#include "ElementALE.h"
#include "Element3D.h"
#include "MathVector.h"


class Element3D4nodedTetrahedron: public virtual Element, public ElementALE, public Element3D
{
  public:

    Element3D4nodedTetrahedron(void);

    virtual ~Element3D4nodedTetrahedron();
	  
    int ndm(void) { return 3; }
    int nen(void) { return 4; }
    
    virtual void plotOutline(bool defFlg = false);
    
    virtual void paint(bool defFlg = false);

    virtual bool forDomainType(int);

    virtual void putLabel(char*, bool defFlg = false);
    
    virtual void contourPlot(int, int, int, double, double, bool defFlg = true);

    virtual double volume(bool init = false);

    virtual int  nFaces(void) { return 4; }

    virtual void giveFace3D(int, Vector<int> &);

    virtual int  calcStiffnessAndResidualMesh(void);

    virtual int  nBasicFacesPerFace(void) { return 1; }

    virtual void defineBasicFace(int, int, int*, unsigned int*);

  private:

	  
};

#endif

