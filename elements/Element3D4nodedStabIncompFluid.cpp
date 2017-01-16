
#include "Element3D4nodedStabIncompFluid.h"
#include "Debug.h"
#include "Plot.h"
#include "ElementGroup.h"
#include "FunctionsElement.h"
#include "PropertyTypeEnum.h"
#include "MpapTime.h"


extern MpapTime   mpapTime;
extern Plot plot;


using namespace std;



Element3D4nodedStabIncompFluid::Element3D4nodedStabIncompFluid(void)
{
  if (debug) cout << " constructor Element3D4nodedStabIncompFluid\n\n";

  return;
}




Element3D4nodedStabIncompFluid::~Element3D4nodedStabIncompFluid()
{
  if (debug) cout << " destructor Element3D4nodedStabIncompFluid\n\n";

  return;
}







bool Element3D4nodedStabIncompFluid::forDomainType(int domType)
{
  switch (domType)
  {
    case MESH:  return true;
		
    case FLUID: return true;
		
    default:    return false;
  }
}





int Element3D4nodedStabIncompFluid::calcStiffnessAndResidual(void)
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
 
  element3dstabincompfluid_(elmDat,timDat,xl,ul,dt,s,p,&ndf_,&ndm_,&nen_,&nelm,&isw);
 
 
  return 0;
}










void Element3D4nodedStabIncompFluid::projectVorticity(int indx)
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
         isw     = 10+indx;

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
  
  element3dstabincompfluid_(elmDat,timDat,xl,ul,dt,s,p,&ndf_,&ndm_,&nen_,&nelm,&isw);
 
 
  return ;
}










void Element3D4nodedStabIncompFluid::projectError(int indx)
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
         isw     = 14;

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
  
  element3dstabincompfluid_(elmDat,timDat,xl,ul,dt,s,p,&ndf_,&ndm_,&nen_,&nelm,&isw);
 
  return;
}



int Element3D4nodedStabIncompFluid::calcMeshDerivatives(void)
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
         isw     = 9;

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
  for (i=0; i<nstu*nstx; i++) s[i] = 0.;
  
  element3dstabincompfluid_(elmDat,timDat,xl,ul,dt,s,p,&ndf_,&ndm_,&nen_,&nelm,&isw);


  //cout << "calcMeshDerivatives\n\n";
  //for (i=0; i<nstu; i++)
  //{
  //  for (j=0; j<nstx; j++) printf(" %12.5g",s[i+j*nstu]); cout << "\n";
  //}
  //cout << "\n";
  //for (j=0; j<nstu; j++) printf(" %12.5g",p[j]); cout << "\n";
  //cout << "\n";

  //exit(0);*/
  
  return 0;
}


void Element3D4nodedStabIncompFluid::diffMeshDerivTest(double ddd, int dig, int dig2, bool gfrmt)
{
  int    i, jnd, jj, j, k, 
         nstu = nen()*ndf(), 
         nstx = nen()*ndm(), 
         niv = nGaussPoints() * nivGP();
	
  ElementGroup *eG = (ElementGroup*) elemGrp;
	
  double *s     = eG->dom->s,
         *p     = eG->dom->p,
         dd[6]  = {-3.*ddd, -2.*ddd, -ddd, +ddd, +2.*ddd, +3.*ddd },
	 *r     = new double[6*nstu],
	 *sdiff = new double[nstu*nstx],
	 smax, sdmax;

  COUT << "element type: " << ((ElementGroup*)elemGrp)->elemProp[ELEMENTTYPE]->name << "\n\n";
 
  for (i=0; i<nivGP()*nGaussPoints(); i++) intVar1[i] = intVar2[i];
 
  // loop over columns of s

  for (jnd=0; jnd<nen(); jnd++) // nodes
  {
    for (jj=0; jj<ndm(); jj++)  // dof
    {
      j = jnd*ndm() + jj;
	    
      // loop over perturbations
	    
      for (k=0; k<6; k++)
      {
        // apply pertubation
	      
        x(jnd+1,jj+1) += dd[k];

        v(jnd+1,jj+1) += dd[k] * eG->dom->td[20];

        eG->dom->updateIterStep();
	
        // calculate residual

	calcStiffnessAndResidual();
	//calcMeshDerivatives();

        for (i=0; i<niv; i++) intVar2[i] = intVar1[i];
  
        // remove pertubation
	
        x(jnd+1,jj+1) -= dd[k];
	
        v(jnd+1,jj+1) -= dd[k] * eG->dom->td[20];

        eG->dom->updateIterStep();
	
        // loop over rows of s, store residual

        for (i=0; i<nstu; i++) r[i*6+k] = p[i];
      }

      // loop over rows of s

      for (i=0; i<nstu; i++)
		
        sdiff[j*nstu+i] = ( +       r[i*6+0]
                            -  9. * r[i*6+1]
                            + 45. * r[i*6+2]
                            - 45. * r[i*6+3]
                            +  9. * r[i*6+4]
                            -       r[i*6+5] ) / (60. * ddd);
    }
  }
 
  // calculate stiffness

  calcMeshDerivatives();

  for (i=0; i<niv; i++) intVar2[i] = intVar1[i];
  
  // compare the matrices

  prgCompareTwoSimpleMatrices(sdiff,                         // matrix 1
		              s,                             // matrix 2
		              "numerical differentiation",   // title matrix 1 
			      "analytical calculation",      // title matrix 2
			      "numerical - analytical",      // title matrix 1 - 2
			      nstu, nstx,                    // matrix dimension
			      dig,dig2,gfrmt,                // format
			      0,                             // indentation
			      false,                         // interactive
			      false);                        // row/column numbers
  delete [] r;
  delete [] sdiff;

  return;
}











