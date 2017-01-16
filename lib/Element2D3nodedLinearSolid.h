
#ifndef incl_Element2D3nodedLinearSolid_h
#define incl_Element2D3nodedLinearSolid_h


#include "ElementSolid.h"
#include "Element2D3nodedTriangle.h"



class Element2D3nodedLinearSolid: public Element2D3nodedTriangle, public ElementSolid
{
  public:

    Element2D3nodedLinearSolid(void);

    virtual ~Element2D3nodedLinearSolid();
	  
    virtual bool forDomainType(int);

    virtual int stressStrainState(void);

    virtual int finiteStrain(void);

    virtual int nGaussPoints(void){return 1;}
    
    int calcStiffnessAndResidual(void);

    virtual int ndf(void) { return 2; }

    virtual void projectIntVar(int);
    
    virtual void projectStress(int);

    virtual void plotGaussPoints(int, bool defFlg = false);

  private:

};

#endif

