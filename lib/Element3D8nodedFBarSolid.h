
#ifndef incl_Element3D8nodedFBarSolid_h
#define incl_Element3D8nodedFBarSolid_h


#include "ElementSolid.h"
#include "Element3D8nodedBrick.h"



class Element3D8nodedFBarSolid: public Element3D8nodedBrick, public ElementSolid
{
  public:

    Element3D8nodedFBarSolid(void);

    virtual ~Element3D8nodedFBarSolid();
	  
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

