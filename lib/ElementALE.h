
#ifndef incl_ElementALE_h
#define incl_ElementALE_h


#include "Element.h"


class ElementALE: public virtual Element
{
  public:

    //ElementALE(void);

    //virtual ~ElementALE();

    virtual int ndfm(int what) { if (what == MSH) return ndm(); else return ndf(); }
    
    virtual int calcStiffnessAndResidualMesh(void) = 0;
    
    virtual void diffAleTest(double,int,int,bool);
	    
  private:

	  
};

#endif

