
#include "Element2D9nodedQuadraticSolid.h"
#include "Plot.h"
#include "ElementGroup.h"
#include "FunctionsElement.h"
#include "FunctionsProgram.h"
#include "PropertyTypeEnum.h"
#include "Mesh.h"
#include "Solid.h"


extern Plot plot;


using namespace std;



Element2D9nodedQuadraticSolid::Element2D9nodedQuadraticSolid(void)
{
  return;
}




Element2D9nodedQuadraticSolid::~Element2D9nodedQuadraticSolid()
{
  return;
}







bool Element2D9nodedQuadraticSolid::forDomainType(int domType)
{
  switch (domType)
  {
    case MESH:          return true;
		
    case SOLID:         return true;
		
    case MICROCELLWULF: return true;
		
    case MICROCELL:     return true;
		
    default:            return false;
  }
}





int Element2D9nodedQuadraticSolid::stressStrainState(void)
{
  if (((ElementGroup*)elemGrp)->elemProp[ELEMENTTYPE]->data.n < 3)

    prgError(1,"Element2D9nodedQuadraticSolid::stressStrainState","element type data empty!");
	
  return roundToInt(((ElementGroup*)elemGrp)->elemProp[ELEMENTTYPE]->data[2]);
}





int Element2D9nodedQuadraticSolid::finiteStrain(void)
{
  if (((ElementGroup*)elemGrp)->elemProp[ELEMENTTYPE]->data.n < 2)

    prgError(1,"Element2D9nodedQuadraticSolid::finiteStrain","element type data empty!");
	
  return roundToInt( ((ElementGroup*)elemGrp)->elemProp[ELEMENTTYPE]->data[1]);
}





int Element2D9nodedQuadraticSolid::nGaussPoints(void)
{
  if (((ElementGroup*)elemGrp)->elemProp[ELEMENTTYPE]->data.n < 1)

    prgError(1,"Element2D9nodedQuadraticSolid::nGaussPoints","element type data empty!");
	
  int ngp = roundToInt(((ElementGroup*)elemGrp)->elemProp[ELEMENTTYPE]->data[0]);

  if (ngp!=1 && ngp!=4 && ngp!=9) 
	  
    prgError(2,"Element2D9nodedQuadraticSolid::nGaussPoints","invalid number of Gauss points!");

  return ngp;
}





int Element2D9nodedQuadraticSolid::calcStiffnessAndResidual(void)
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
    else prgError(1,"Element2D9nodedQuadraticSolid::calcStiffnessAndResidual",
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

 /* for (i=0; i<nstu; i++)
  {
    for (j=0; j<nstu; j++) printf(" %12.5g",s[i+j*nstu]); cout << "\n";
  }
  cout << "\n";
  for (j=0; j<nstu; j++) printf(" %12.5g",p[j]); cout << "\n";
  cout << "\n";*/
  
  return errCode;
}





void Element2D9nodedQuadraticSolid::projectIntVar(int indx)
{
  projectToNodes(intVar2, indx, ((ElementGroup*)elemGrp)->nivGP);

  return;
}






void Element2D9nodedQuadraticSolid::projectStress(int indx)
{
  ElementGroup *eG = (ElementGroup*) elemGrp;
	
  double *s      = eG->dom->s,
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
    { prgWarning(1,"Element2D9nodedQuadraticSolid::projectStress","nstu*nstu < ngp * 4!"); return; }
  
  element2dsolid_(elmDat,matDat,timDat,xl,ul,intVar1,intVar2,s,s,
                  &ndf_,&ndm_,&nen_,&nelm,&nmat,&nivGP,&matId,&isw,
                  &errCode,(Element*) this);
 
  projectToNodes(s, indx, 4);

  return;
}






void Element2D9nodedQuadraticSolid::projectToNodes(double *dat, int indx, int nn)
{
  ElementGroup *eG = (ElementGroup*) elemGrp;
  
  int i, l, ll[10], ii = indx, 
      ngp = roundToInt(eG->elemProp[ELEMENTTYPE]->data[0]);

  double fct4[4] = { 1.86603e0, -0.5e0, 0.133975e0, -0.5e0 },
         *p = eG->dom->p;
	
  if (nn < 1) 
    { COUT << "nn = 1!  Element2D9nodedQuadraticSolid::projectToNodes aborted!\n\n"; return; }
  
  if (nn < ii) ii = nn; if (ii < 1) ii = 1;   ii--;
 
  for (i=0; i<nen(); i++) p[i] = 0.;
  
  switch (ngp)
  {

    case  9: cout << "    What shall we do with the drunken sailor? \n\n";
	    
	     break;

    default: prgError(1,"Element2D9nodedQuadraticSolid::projectToNodes","invalid ngp");
  }
	  
  return;
}






void Element2D9nodedQuadraticSolid::plotGaussPoints(int num, bool defFlg)
{
  double d = (plot.dAct[0] + plot.dAct[1]) * .0055,

         X[18], xi[18], wgp[9], shp[27], detJ, xgp[3];
 
  int    l, i, ngp = nGaussPoints(), nen_ = 8;

  if (ngp > 9) prgError(1,"Element2D9nodedQuadraticSolid::plotGaussPoints","ngp > 9 !");

  if (ngp != GpDatId.n) return;
  
  if (defFlg) for (i=0; i<9; i++) { X[i+i] = x (i+1,1); X[i+i+1] = x (i+1,2); }
  else        for (i=0; i<9; i++) { X[i+i] = x0(i+1,1); X[i+i+1] = x0(i+1,2); }
  
  compxigp2d_(xi,wgp,&ngp);

  for (l=0; l<ngp; l++) 
  {
    compshp2d_(shp,&detJ,&(xi[l+l]),X,&nen_);  // we only use the shape function values

    xgp[0] = 0.; for (i=0; i<9; i++) xgp[0] += shp[i*3+2] * X[i+i];	    
    xgp[1] = 0.; for (i=0; i<9; i++) xgp[1] += shp[i*3+2] * X[i+i+1];	    
    
    //cout << " " << xgp[0] << ", " << xgp[1] << "\n"; cout << "\n";

    if (num == 1)  plot.point(xgp,d,GpDatId[l]+1); else plot.point(xgp,d);
  }

  return;
}










void Element2D9nodedQuadraticSolid::projectError(int)
{
  ElementGroup *eG = (ElementGroup*) elemGrp;
	
  double *s      = eG->dom->s,
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
         isw     = 8,
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
    { prgWarning(1,"Element2D9nodedQuadraticSolid::projectError","nstu*nstu < ngp!"); return; }
  
  element2dsolid_(elmDat,matDat,timDat,xl,ul,intVar1,intVar2,s,s,
                  &ndf_,&ndm_,&nen_,&nelm,&nmat,&nivGP,&matId,&isw,
                  &errCode,(Element*) this);
 
  projectToNodes(s, 1, 1);

  return;
}

