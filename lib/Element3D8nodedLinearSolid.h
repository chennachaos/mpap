
#ifndef incl_Element3D8nodedLinearSolid_h
#define incl_Element3D8nodedLinearSolid_h


#include "ElementSolid.h"
#include "Element3D8nodedBrick.h"



class Element3D8nodedLinearSolid: public Element3D8nodedBrick, public ElementSolid
{
  public:

    Element3D8nodedLinearSolid(void);

    virtual ~Element3D8nodedLinearSolid();
	  
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

