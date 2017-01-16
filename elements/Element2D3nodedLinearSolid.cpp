
#include "Element2D3nodedLinearSolid.h"
#include "Plot.h"
#include "ElementGroup.h"
#include "FunctionsElement.h"
#include "FunctionsProgram.h"
#include "PropertyTypeEnum.h"
#include "Mesh.h"
#include "Solid.h"


extern Plot plot;


using namespace std;



Element2D3nodedLinearSolid::Element2D3nodedLinearSolid(void)
{
  return;
}




Element2D3nodedLinearSolid::~Element2D3nodedLinearSolid()
{
  return;
}







bool Element2D3nodedLinearSolid::forDomainType(int domType)
{
  switch (domType)
  {
    case MESH:          return true;
		
    case SOLID:         return true;
		
    case MICROCELL: return true;

    case MICROCELLWULF: return true;
		
    default:            return false;
  }
}





int Element2D3nodedLinearSolid::stressStrainState(void)
{
  if (((ElementGroup*)elemGrp)->elemProp[ELEMENTTYPE]->data.n < 3)

    prgError(1,"Element2D3nodedLinearSolid::stressStrainState","element type data empty!");
	
  return roundToInt(((ElementGroup*)elemGrp)->elemProp[ELEMENTTYPE]->data[2]);
}





int Element2D3nodedLinearSolid::finiteStrain(void)
{
  if (((ElementGroup*)elemGrp)->elemProp[ELEMENTTYPE]->data.n < 2)

    prgError(1,"Element2D3nodedLinearSolid::finiteStrain","element type data empty!");
	
  return roundToInt( ((ElementGroup*)elemGrp)->elemProp[ELEMENTTYPE]->data[1]);
}





int Element2D3nodedLinearSolid::calcStiffnessAndResidual(void)
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
         isw     = 4,
	 errCode = 0;

  if (eG->dom->tis != 0) 
  {
    if (((Solid*)eG->dom)->lumpType == 0) isw = 6;
    else prgError(1,"Element2D3nodedLinearSolid::calcStiffnessAndResidual",
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

  
  // zero element residual and stiffness
  
  for (i=0; i<nstu;      i++) p[i] = 0.;
  for (i=0; i<nstu*nstu; i++) s[i] = 0.;
  
  element2dsolid_(elmDat,matDat,timDat,xl,ul,intVar1,intVar2,s,p,
                           &ndf_,&ndm_,&nen_,&nelm,&nmat,&nivGP,&matId,&isw,
	                   &errCode,(Element*) this);
 
//  for (i=0; i<nstu; i++)
//  {
//    for (j=0; j<nstu; j++) printf(" %12.5g",s[i+j*nstu]); cout << "\n";
//  }
//  cout << "\n";
//  for (j=0; j<nstu; j++) printf(" %12.5g",p[j]); cout << "\n";
//  cout << "\n";
  
  return errCode;
}








void Element2D3nodedLinearSolid::projectIntVar(int indx)
{
  ElementGroup *eG = (ElementGroup*) elemGrp;

  if (eG->nivGP == 0) 
  {
    prgWarning(1,"Element2D3nodedLinearSolid::projectIntVar","no internal variables!");
    return;
  }
  
  int i = indx;
  
  double *p = eG->dom->p;
 
  if (i > eG->nivGP) i = eG->nivGP;
  if (i < 1)         i = 1;
	 
  p[0] = intVar2[--i];

  p[1] = p[0];
  p[2] = p[0];
  
  return;
}






void Element2D3nodedLinearSolid::projectStress(int indx)
{
  ElementGroup *eG = (ElementGroup*) elemGrp;
	
  double *p      = eG->dom->p,
         *s      = eG->dom->s,
         *xl     = eG->dom->xl,
         *ul     = eG->dom->ul,
         *elmDat = &(eG->elemProp[ELEMENTTYPE]->data[0]),
         *matDat = &(eG->elemProp[MATERIAL]->data[0]),
         *timDat = eG->dom->td;
 
  int    i, j, nivEl,
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
	 ngp     = roundToInt(elmDat[0]),
         nelm    = eG->elemProp[ELEMENTTYPE]->data.n,
         nmat    = eG->elemProp[MATERIAL]->data.n,
	 matId   = eG->elemProp[MATERIAL]->id + 1,
     isw     = 7,
	 errCode = 0;
  
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

  for (i=0; i<nstu*nstu; i++) s[i] = 0.;

  nivEl = ngp * nivGP;
  
  for (i=0; i<nivEl; i++) intVar2[i] = intVar1[i];

  if (nstu*nstu < ngp * 4) 
    { prgWarning(1,"Element2D3nodedLinearSolid::projectStress","nstu*nstu < ngp * 4!"); return; }
  
  element2dsolid_(elmDat,matDat,timDat,xl,ul,intVar1,intVar2,s,s,
		             &ndf_,&ndm_,&nen_,&nelm,&nmat,&nivGP,&matId,&isw,
			     &errCode,(Element*) this);
  i = --indx; 
  if (i<0) i = 0;
  if (i>3) i = 3;
  
  p[0] = s[i];
  p[1] = p[0];
  p[2] = p[1];
  
  return ;
}






void Element2D3nodedLinearSolid::plotGaussPoints(int num, bool defFlg)
{
  double d = (plot.dAct[0] + plot.dAct[1]) * .0055, X[2] = { 0., 0.};
 
  int    i, ngp = 1;

  if (ngp != GpDatId.n) return;
  
  if (defFlg) for (i=0; i<3; i++) { X[0] += x (i+1,1); X[1] += x (i+1,2); }
  else        for (i=0; i<3; i++) { X[0] += x0(i+1,1); X[1] += x0(i+1,2); }
  
  X[0] /= 3.;
  X[1] /= 3.;
  
  if (num == 1)  plot.point(X,d,GpDatId[0]+1); else plot.point(X,d);

  return;
}





