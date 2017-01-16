
#ifndef incl_Element3D10nodedQuadraticSolid_h
#define incl_Element3D10nodedQuadraticSolid_h


#include "ElementSolid.h"
#include "Element3D10nodedTetrahedron.h"



class Element3D10nodedQuadraticSolid: public Element3D10nodedTetrahedron, public ElementSolid
{
  public:

    Element3D10nodedQuadraticSolid(void);

    virtual ~Element3D10nodedQuadraticSolid();
	  
    virtual bool forDomainType(int);

    virtual int finiteStrain(void);

	virtual int nGaussPoints(void){return 4;}
    
    int calcStiffnessAndResidual(void);

    virtual int ndf(void) { return 3; }

    virtual void projectIntVar(int);
    
    virtual void projectStress(int);

    virtual void plotGaussPoints(int, bool defFlg = false);

  private:

    void projectToNodes(double*, int, int);
    
};

#endif

