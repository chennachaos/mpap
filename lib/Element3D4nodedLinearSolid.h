
#ifndef incl_Element3D4nodedLinearSolid_h
#define incl_Element3D4nodedLinearSolid_h


#include "ElementSolid.h"
#include "Element3D4nodedTetrahedron.h"



class Element3D4nodedLinearSolid: public Element3D4nodedTetrahedron, public ElementSolid
{
  public:

    Element3D4nodedLinearSolid(void);

    virtual ~Element3D4nodedLinearSolid();
	  
    virtual bool forDomainType(int);

    virtual int finiteStrain(void);

    virtual int nGaussPoints(void){return 1;}
    
    int calcStiffnessAndResidual(void);

    virtual int ndf(void) { return 3; }

    virtual void projectIntVar(int);
    
    virtual void projectStress(int);

    virtual void plotGaussPoints(int, bool defFlg = false);

  private:

    void projectToNodes(double*, int, int);
    
};

#endif

