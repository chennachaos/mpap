
#ifndef incl_ElementSolid_h
#define incl_ElementSolid_h


#include "Element.h"


class ElementSolid: public virtual Element
{
  public:

    ElementSolid(void);

    virtual ~ElementSolid();

    virtual int finiteStrain(void);
    
    virtual int stressStrainState(void);
    
    virtual int nGaussPoints(void);
    
    virtual int nivGP(void);

    virtual void initialiseIntVar(void);

    virtual void setGaussPointDataId(void);

  private:

	  
};

#endif

