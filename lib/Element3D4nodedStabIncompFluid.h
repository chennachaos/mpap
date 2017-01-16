#ifndef incl_Element3D4nodedStabIncompFluid_h
#define incl_Element3D4nodedStabIncompFluid_h


#include "Element3D4nodedTetrahedron.h"



class Element3D4nodedStabIncompFluid: public Element3D4nodedTetrahedron
{
  public:

    Element3D4nodedStabIncompFluid(void);

    virtual ~Element3D4nodedStabIncompFluid();
	  
    virtual bool forDomainType(int);

    int calcStiffnessAndResidual(void);

    void projectVorticity(int);

    virtual int ndf(void) { return 4; }
    
    virtual void projectError(int);

    virtual int calcMeshDerivatives(void);

    virtual void diffMeshDerivTest(double,int,int,bool);

  private:

	  
};

#endif



