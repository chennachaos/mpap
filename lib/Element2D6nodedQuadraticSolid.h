
#ifndef incl_Element2D6nodedQuadraticSolid_h
#define incl_Element2D6nodedQuadraticSolid_h


#include "ElementSolid.h"
#include "Element2D6nodedTriangle.h"



class Element2D6nodedQuadraticSolid: public Element2D6nodedTriangle, public ElementSolid
{
  public:

    Element2D6nodedQuadraticSolid(void);

    virtual ~Element2D6nodedQuadraticSolid();
	  
    virtual bool forDomainType(int);

    virtual int stressStrainState(void);

    virtual int finiteStrain(void);

    virtual int nGaussPoints(void){return 3;}
    
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

