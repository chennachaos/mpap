
#ifndef incl_Element2D3nodedLine_h
#define incl_Element2D3nodedLine_h


#include "Element.h"
#include "Element2D.h"



class Element2D3nodedLine: public virtual Element, public Element2D
{
  public:

    Element2D3nodedLine(void);

    virtual ~Element2D3nodedLine();
	  
    int ndm(void) { return 2; }
    int nen(void) { return 3; }
    
    virtual void plotOutline(bool defFlg = false);
    
    virtual void paint(bool defFlg = false);

    virtual bool forDomainType(int);

    virtual void putLabel(char*, bool defFlg = false);

    virtual void givePlotSequence2D(Vector<int> &);

    virtual bool isClosed(void) { return false; }
    
  private:

	  
};

#endif

