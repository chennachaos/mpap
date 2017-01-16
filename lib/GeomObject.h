
#ifndef incl_GeomObject_h
#define incl_GeomObject_h


#include "List.h"
#include "GeomTypeEnum.h"


class GeomObject: public ListItem
{ 
  public:

    GeomObject(void) { return; }

    virtual ~GeomObject() { return; }
  
    int dat1;

    virtual bool isType(GeomTypeEnum) = 0;

    Vector<int> nodeListTmp;
  
  private:

};

#endif










