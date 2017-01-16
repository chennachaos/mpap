
#ifndef incl_Element2D8nodedQuadraticSolid_h
#define incl_Element2D8nodedQuadraticSolid_h


#include "ElementSolid.h"
#include "Element2D8nodedQuadrilateral.h"



class Element2D8nodedQuadraticSolid: public Element2D8nodedQuadrilateral, public ElementSolid
{
  public:

    Element2D8nodedQuadraticSolid(void);

    virtual ~Element2D8nodedQuadraticSolid();
	  
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

