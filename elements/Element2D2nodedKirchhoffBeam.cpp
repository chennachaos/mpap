
#include "Element2D2nodedKirchhoffBeam.h"
#include "Plot.h"
#include "ElementGroup.h"
#include "FunctionsElement.h"
#include "FunctionsProgram.h"
#include "PropertyTypeEnum.h"
#include "Mesh.h"
#include "Solid.h"


extern Plot plot;


using namespace std;



Element2D2nodedKirchhoffBeam::Element2D2nodedKirchhoffBeam(void)
{
  return;
}




Element2D2nodedKirchhoffBeam::~Element2D2nodedKirchhoffBeam()
{
  return;
}







bool Element2D2nodedKirchhoffBeam::forDomainType(int domType)
{
  switch (domType)
  {
    case MESH:      return true;
		
    case SOLID:     return true;
		
    default:        return false;
  }
}







int Element2D2nodedKirchhoffBeam::calcStiffnessAndResidual(void)
{
  ElementGroup *eG = (ElementGroup*) elemGrp;
	
  double *s      = eG->dom->s,
         *p      = eG->dom->p,
         *xl     = eG->dom->xl,
         *ul     = eG->dom->ul,
         *elmDat = &(eG->elemProp[ELEMENTTYPE]->data[0]),
         *matDat = &(eG->elemProp[MATERIAL]->data[0]),
         *timDat = eG->dom->td;
 
  int    i, j,
	 nstu    = ndf()*nen(),
	 nstu2   = nstu  + nstu,
	 nstu3   = nstu2 + nstu,
	 nstu4   = nstu3 + nstu,
	 nstu5   = nstu4 + nstu,
	 nstx    = ndm()*nen(),
	 nstx2   = nstx  + nstx,
	 ndf_    = ndf(),
         ndm_    = ndm(),
         nen_    = nen(),
         nivGP   = eG->nivGP,
         nelm    = eG->elemProp[ELEMENTTYPE]->data.n,
         nmat    = eG->elemProp[MATERIAL]->data.n,
	 matId   = eG->elemProp[MATERIAL]->id + 1,
         isw     = 4;

  if (eG->dom->tis != 0) 
  {
    if (((Solid*)eG->dom)->lumpType == 0) isw = 6;
    else prgError(1,"Element2D2nodedKirchhoffBeam::calcStiffnessAndResidual",
		    "mass lumping not available yet!");
  }

  

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
      xl[i*ndm_+j+nstx2] = x0(i+1,j+1);
    }
  }

  for (i=0; i<nstu;      i++) p[i] = 0.;
  for (i=0; i<nstu*nstu; i++) s[i] = 0.;
 
  element2dkirchhoffbeam_(elmDat,matDat,timDat,xl,ul,intVar1,intVar2,s,p,
		                  &ndf_,&ndm_,&nen_,&nelm,&nmat,&nivGP,&matId,&isw,
				  (Element*) this);
 
//  for (i=0; i<nstu; i++)
//  {
//    for (j=0; j<nstu; j++) printf(" %12.5g",s[i+j*nstu]); cout << "\n";
//  }
//  cout << "\n";
//  for (j=0; j<nstu; j++) printf(" %12.5g",p[j]); cout << "\n";
//  cout << "\n";
  
  return 0;
}

