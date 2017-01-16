
#ifndef incl_Element2D2nodedGeomExSmallStrainBeam_h
#define incl_Element2D2nodedGeomExSmallStrainBeam_h


#include "ElementSolid.h"
#include "Element2D2nodedLine.h"



class Element2D2nodedGeomExSmallStrainBeam: public Element2D2nodedLine, public ElementSolid
{
  public:

    Element2D2nodedGeomExSmallStrainBeam(void);

    virtual ~Element2D2nodedGeomExSmallStrainBeam();
	  
    virtual bool forDomainType(int);

    int calcStiffnessAndResidual(void);

    virtual int nGaussPoints(void) { return 1; }
    
    virtual int ndf(void) { return 3; }

    virtual double volume(bool init = false) { return 0.; }

    virtual void contourPlot(int, int, int, double, double, bool defFlg = true) { return; }
   
    virtual void projectIntVar(int) { return; }
    
    virtual void projectStress(int) { return; }
    
  private:

	  
};

#endif

