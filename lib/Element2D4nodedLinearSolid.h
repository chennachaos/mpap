
#ifndef incl_Element2D4nodedLinearSolid_h
#define incl_Element2D4nodedLinearSolid_h


#include "ElementSolid.h"
#include "Element2D4nodedQuadrilateral.h"



class Element2D4nodedLinearSolid: public Element2D4nodedQuadrilateral, public ElementSolid
{
  public:

    Element2D4nodedLinearSolid(void);

    virtual ~Element2D4nodedLinearSolid();
	  
    virtual bool forDomainType(int);

    virtual int stressStrainState(void);

    virtual int finiteStrain(void);

    virtual int nGaussPoints(void);
    
    int calcStiffnessAndResidual(void);

    virtual int ndf(void) { return 2; }

    virtual void projectIntVar(int);
    
    virtual void projectStress(int);

    virtual void projectError(int);

    virtual void plotGaussPoints(int, bool defFlg = false);

  private:

    void projectToNodes(double*, int, int);
    
};

#endif

