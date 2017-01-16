
#include "Element1D2nodedLinearSolidALE.h"
#include "Plot.h"
#include "ElementGroup.h"
#include "FunctionsElement.h"
#include "FunctionsProgram.h"
#include "PropertyTypeEnum.h"
#include "Mesh.h"
#include "MpapTime.h"


extern Plot     plot;
extern MpapTime mpapTime;


using namespace std;



Element1D2nodedLinearSolidALE::Element1D2nodedLinearSolidALE(void)
{
  return;
}




Element1D2nodedLinearSolidALE::~Element1D2nodedLinearSolidALE()
{
  return;
}







bool Element1D2nodedLinearSolidALE::forDomainType(int domType)
{
  switch (domType)
  {
    case FLUID: return true;
		
    default:    return false;
  }
}







int Element1D2nodedLinearSolidALE::calcStiffnessAndResidual(void)
{
  ElementGroup *eG = (ElementGroup*) elemGrp;
	
  double *s      = eG->dom->s,
         *p      = eG->dom->p,
         *xl     = eG->dom->xl,
         *ul     = eG->dom->ul,
         *elmDat = &(eG->elemProp[ELEMENTTYPE]->data[0]),
         *timDat = eG->dom->td;
 
  int    i, j,
	 nstu    = ndf()*nen(),
	 nstu2   = nstu  + nstu,
	 nstu3   = nstu2 + nstu,
	 nstx    = ndm()*nen(),
	 nstx2   = nstx  + nstx,
	 ndf_    = ndf(),
         ndm_    = ndm(),
         nen_    = nen(),
         nivGP   = eG->nivGP,
         nelm    = eG->elemProp[ELEMENTTYPE]->data.n,
         isw     = 4,
	 errCode = 0;
 
  for (i=0; i<nen_; i++)
  {
    for (j=0; j<ndf_; j++)
    {
      ul[i*ndf_+j]       =  u(i+1,j+1);
      ul[i*ndf_+j+nstu]  = un(i+1,j+1);
      ul[i*ndf_+j+nstu2] = u3(i+1,j+1);
      ul[i*ndf_+j+nstu3] = u4(i+1,j+1);
    }
    for (j=0; j<ndm_; j++) 
    {
      xl[i*ndm_+j]       =  x(i+1,j+1);
      xl[i*ndm_+j+nstx]  = xn(i+1,j+1);
      xl[i*ndm_+j+nstx2] = x0(i+1,j+1);
    }
  }
  
  // zero element residual and stiffness
  
  for (i=0; i<nstu;      i++) p[i] = 0.;
  for (i=0; i<nstu*nstu; i++) s[i] = 0.;
  
  element1dsolidale_(elmDat,timDat,xl,ul,s,p);

  return 0;
 
  for (i=0; i<nstu; i++)
  {
    for (j=0; j<nstu; j++) printf(" %12.5g",s[i+j*nstu]); cout << "\n";
  }
  cout << "\n";
  for (j=0; j<nstu; j++) printf(" %12.5g",p[j]); cout << "\n";
  cout << "\n";

  return 0;
}

