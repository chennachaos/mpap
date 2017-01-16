
#ifndef incl_ElementGroup_h
#define incl_ElementGroup_h


#include "MyString.h"
#include "List.h"
#include "Element.h"
#include "MathVector.h"
#include "Domain.h"
#include "PropertyItem.h"



class ElementGroup: public ListItem
{
  public:

    ElementGroup(void);

    ElementGroup(Domain*);

    ~ElementGroup();
	  
    List<Element> elem;

    Domain *dom;
    
    int ndf, ndm, nen, nivGP, nGP;

    bool closed;
    
    int nElemProp;
    
    PropertyItem **elemProp;

    Vector<int> geomObj;

    void addElemProp(int, PropertyItem*);
    
  private:
    
};

#endif

