
#ifndef incl_Element2D3nodedTriangle_h
#define incl_Element2D3nodedTriangle_h


#include "Element.h"
#include "ElementALE.h"
#include "Element2D.h"
#include "MathVector.h"



class Element2D3nodedTriangle: public virtual Element, public ElementALE, public Element2D
{
  public:

    Element2D3nodedTriangle(void);

    virtual ~Element2D3nodedTriangle();
	  
    virtual int ndm(void) { return 2; }
    virtual int nen(void) { return 3; }
    
    virtual void plotOutline(bool defFlg = false);
    
    virtual void paint(bool defFlg = false);

    virtual bool forDomainType(int);

    virtual void putLabel(char*, bool defFlg = false);

    virtual void contourPlot(int, int, int, double, double, bool defFlg = true);

    virtual void givePlotSequence2D(Vector<int> &);

    virtual double volume(bool init = false);
    
    virtual double diameter(bool init = false);
    
    virtual void projectVorticity(int);

    virtual void projectGradient(int,int);

    virtual void projectNormOfGradient(int);

    virtual void projectNormOfGradientSquared(int);

    virtual void getGradient(int, double *);

    virtual int  calcStiffnessAndResidualMesh(void);

    virtual bool containsPoint(double*, double*);

  private:

	  
};

#endif

