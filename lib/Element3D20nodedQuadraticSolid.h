
#ifndef incl_Element3D20nodedQuadraticSolid_h
#define incl_Element3D20nodedQuadraticSolid_h


#include "ElementSolid.h"
#include "Element3D20nodedBrick.h"



class Element3D20nodedQuadraticSolid: public Element3D20nodedBrick, public ElementSolid
{
  public:

    Element3D20nodedQuadraticSolid(void);

    virtual ~Element3D20nodedQuadraticSolid();
	  
    virtual bool forDomainType(int);

    virtual int finiteStrain(void);

    virtual int nGaussPoints(void);
    
    int calcStiffnessAndResidual(void);

    virtual int ndf(void) { return 3; }

    virtual void projectIntVar(int);
    
    virtual void projectStress(int);

    virtual void plotGaussPoints(int, bool defFlg = false);

  private:

    void projectToNodes(double*, int, int);
    
};

#endif

