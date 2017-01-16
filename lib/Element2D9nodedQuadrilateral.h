
#ifndef incl_Element2D9nodedQuadrilateral_h
#define incl_Element2D9nodedQuadrilateral_h


#include "Element.h"
#include "Element2D.h"



class Element2D9nodedQuadrilateral: public virtual Element, public Element2D
{
  public:

    Element2D9nodedQuadrilateral(void);

    virtual ~Element2D9nodedQuadrilateral();
	  
    int ndm(void) { return 2; }
    int nen(void) { return 9; }
    
    virtual void plotOutline(bool defFlg = false);
    
    virtual void paint(bool defFlg = false);

    virtual bool forDomainType(int);

    virtual void putLabel(char*, bool defFlg = false);

    virtual void contourPlot(int, int, int, double, double, bool defFlg = true);

    virtual double volume(bool init = false);

    virtual void givePlotSequence2D(Vector<int> &);
    
    virtual bool containsPoint(double*, double*);

    virtual double diameter(bool init = false);
    
  private:

	  
};

#endif

