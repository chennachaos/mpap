
#include "Element2D3nodedLinearPoisson.h"
#include "ElementGroup.h"
#include "FunctionsElement.h"
#include "FunctionsProgram.h"
#include "PropertyTypeEnum.h"


using namespace std;



Element2D3nodedLinearPoisson::Element2D3nodedLinearPoisson(void)
{
  return;
}




Element2D3nodedLinearPoisson::~Element2D3nodedLinearPoisson()
{
  return;
}







bool Element2D3nodedLinearPoisson::forDomainType(int domType)
{
  switch (domType)
  {
    case MESH:               return true;
		
    case FINITEELEMENTBVP:   return true;
		
    case FINITEELEMENTBVPWI: return true;

    default:                 return false;
  }
}


//#include "MathGeom.h"


int Element2D3nodedLinearPoisson::calcStiffnessAndResidual(void)
{
  ElementGroup *eG = (ElementGroup*) elemGrp;
	
  double *s      = eG->dom->s,
         *p      = eG->dom->p,
         *xl     = eG->dom->xl,
         *ul     = eG->dom->ul,
         *elmDat = &(eG->elemProp[ELEMENTTYPE]->data[0]),
         xi[3] = { .333333333333333, .333333333333333, .3333333333333333 },
         shp[9], area, f = elmDat[0];
 
  int    i, j = 3,
         nelm = eG->elemProp[ELEMENTTYPE]->data.n;
 
  for (i=0; i<3; i++)
  {
    ul[i]     = u(i+1,1);
    xl[i+i]   = x(i+1,1);
    xl[i+i+1] = x(i+1,2);
  }

  compshp2d_(shp,&area,xi,xl,&j);

  //i = 2; j = 1;
  //trishp_(xi,xl,&i,&j,&area,shp);

  area *= .5;

  //cout << triangleArea2D(xl,xl+2,xl+4) << " = " << area << "\n";

  for (i=0; i<3; i++)

    for (j=0; j<3; j++)

      s[i+j*3] = (shp[i*3]*shp[j*3] + shp[1+i*3]*shp[1+j*3]) * area;

  for (i=0; i<3; i++)

    p[i] = - s[i+0]*ul[0] - s[i+3]*ul[1] - s[i+6]*ul[2] - f * shp[2+i*3] * area;

  return 0;
}





