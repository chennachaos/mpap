
#ifndef incl_Element2D9nodedQuadraticSolid_h
#define incl_Element2D9nodedQuadraticSolid_h


#include "ElementSolid.h"
#include "Element2D9nodedQuadrilateral.h"



class Element2D9nodedQuadraticSolid: public Element2D9nodedQuadrilateral, public ElementSolid
{
  public:

    Element2D9nodedQuadraticSolid(void);

    virtual ~Element2D9nodedQuadraticSolid();
	  
    virtual bool forDomainType(int);

    virtual int stressStrainState(void);

    virtual int finiteStrain(void);

    virtual int nGaussPoints(void);
    
    int calcStiffnessAndResidual(void);

    virtual int ndf(void) { return 2; }

    virtual void projectIntVar(int);
    
    virtual void projectStress(int);

    virtual void plotGaussPoints(int, bool defFlg = false);

    virtual void projectError(int);

  private:

    void projectToNodes(double*, int, int);

	  
};

#endif

