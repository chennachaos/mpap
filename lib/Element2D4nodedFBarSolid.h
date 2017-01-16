
#ifndef incl_Element2D4nodedFBarSolid_h
#define incl_Element2D4nodedFBarSolid_h


#include "ElementSolid.h"
#include "Element2D4nodedQuadrilateral.h"



class Element2D4nodedFBarSolid: public Element2D4nodedQuadrilateral, public ElementSolid
{
  public:

    Element2D4nodedFBarSolid(void);

    virtual ~Element2D4nodedFBarSolid();
	  
    virtual bool forDomainType(int);

    virtual int stressStrainState(void);

    virtual int finiteStrain(void);

    virtual int nGaussPoints(void);
    
    virtual int calcStiffnessAndResidual(void);

    virtual int ndf(void) { return 2; }

    virtual void projectIntVar(int);
    
    virtual void projectStress(int);

    virtual void plotGaussPoints(int, bool defFlg = false);

  private:

    void projectToNodes(double*, int, int);
    
};

#endif

