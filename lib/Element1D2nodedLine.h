
#ifndef incl_Element1D2nodedLine_h
#define incl_Element1D2nodedLine_h


#include "Element.h"



class Element1D2nodedLine: public virtual Element
{
  public:

    Element1D2nodedLine(void);

    virtual ~Element1D2nodedLine();
	  
    int ndm(void) { return 1; }
    int nen(void) { return 2; }
    
    virtual void plotOutline(bool defFlg = false);
    
    virtual void paint(bool defFlg = false);

    virtual void putLabel(char*, bool defFlg = false);
    
    virtual bool forDomainType(int);

    virtual bool isClosed(void) { return false; }
    
  private:

	  
};

#endif

