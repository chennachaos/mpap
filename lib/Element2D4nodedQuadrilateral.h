
#ifndef incl_Element2D4nodedQuadrilateral_h
#define incl_Element2D4nodedQuadrilateral_h


#include "Element.h"
#include "ElementALE.h"
#include "Element2D.h"
#include "MathVector.h"


class Element2D4nodedQuadrilateral: public virtual Element, public ElementALE, public Element2D
{
  public:

    Element2D4nodedQuadrilateral(void);

    virtual ~Element2D4nodedQuadrilateral();
	  
    int ndm(void) { return 2; }
    int nen(void) { return 4; }
    
    virtual void plotOutline(bool defFlg = false);
    
    virtual void paint(bool defFlg = false);

    virtual bool forDomainType(int);

    virtual void putLabel(char*, bool defFlg = false);
    
    virtual void contourPlot(int, int, int, double, double, bool defFlg = true);

    virtual double volume(bool init = false);

    virtual double diameter(bool init = false);

    virtual void givePlotSequence2D(Vector<int> &);
    
    virtual int  calcStiffnessAndResidualMesh(void);

    virtual bool containsPoint(double*, double*);

    virtual void getDistLoadFact(Vector<double> &, Vector<int> &);

  private:
	  
};

#endif

