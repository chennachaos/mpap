
#include "Element2D3nodedStabIncompHighReFluid.h"
#include "Debug.h"
#include "Plot.h"
#include "ElementGroup.h"
#include "FunctionsElement.h"
#include "PropertyTypeEnum.h"
#include "MpapTime.h"


extern MpapTime   mpapTime;
extern Plot plot;


using namespace std;



Element2D3nodedStabIncompHighReFluid::Element2D3nodedStabIncompHighReFluid(void)
{
  if (debug) cout << " constructor Element2D3nodedStabIncompHighReFluid\n\n";

  return;
}




Element2D3nodedStabIncompHighReFluid::~Element2D3nodedStabIncompHighReFluid()
{
  if (debug) cout << " destructor Element2D3nodedStabIncompHighReFluid\n\n";

  return;
}







bool Element2D3nodedStabIncompHighReFluid::forDomainType(int domType)
{
  switch (domType)
  {
    case MESH:  return true;
		
    case FLUID: return true;
		
    default:    return false;
  }
}





int Element2D3nodedStabIncompHighReFluid::calcStiffnessAndResidual(void)
{
  ElementGroup *eG = (ElementGroup*) elemGrp;
	
  double *s      = eG->dom->s,
         *p      = eG->dom->p,
         *xl     = eG->dom->xl,
         *ul     = eG->dom->ul,
         *elmDat = &(eG->elemProp[ELEMENTTYPE]->data[0]),
         *timDat = eG->dom->td,
         *dt     = &(mpapTime.dt);
 
  int    i, j,
	 nstu    = ndf()*nen(),
	 nstu2   = nstu  + nstu,
	 nstu3   = nstu2 + nstu,
	 nstu4   = nstu3 + nstu,
	 nstu5   = nstu4 + nstu,
	 nstx    = ndm()*nen(),
	 ndf_    = ndf(),
	 nstx2   = nstx  + nstx,
	 nstx3   = nstx2 + nstx,
         ndm_    = ndm(),
         nen_    = nen(),
         nelm    = eG->elemProp[ELEMENTTYPE]->data.n,
         isw     = 5;

  for (i=0; i<nen_; i++)
  {
    for (j=0; j<ndf_; j++)
    {
      ul[i*ndf_+j]       =  u(i+1,j+1);
      ul[i*ndf_+j+nstu]  = un(i+1,j+1);
      ul[i*ndf_+j+nstu2] = u3(i+1,j+1);
      ul[i*ndf_+j+nstu3] = u4(i+1,j+1);
      ul[i*ndf_+j+nstu4] = u5(i+1,j+1);
      ul[i*ndf_+j+nstu5] = u6(i+1,j+1);
    }
    for (j=0; j<ndm_; j++) 
    {
      xl[i*ndm_+j]       =  x(i+1,j+1);
      xl[i*ndm_+j+nstx]  = xn(i+1,j+1);
      xl[i*ndm_+j+nstx2] =  v(i+1,j+1);
      xl[i*ndm_+j+nstx3] = vn(i+1,j+1);
    }
  }

  
  // zero element residual and stiffness
  
  for (i=0; i<nstu;      i++) p[i] = 0.;
  for (i=0; i<nstu*nstu; i++) s[i] = 0.;
  
  element2dstabincomphighrefluid_(elmDat,timDat,xl,ul,dt,s,p,&ndf_,&ndm_,&nen_,&nelm,&isw);
 
/*  for (i=0; i<nstu; i++)
  {
    for (j=0; j<nstu; j++) printf(" %12.5g",s[i+j*nstu]); cout << "\n";
  }
  cout << "\n";
  for (j=0; j<nstu; j++) printf(" %12.5g",p[j]); cout << "\n";
  cout << "\n";

  exit(0);*/
  
  return 0;
}

