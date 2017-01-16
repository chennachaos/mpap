
#include "Element2D3nodedPressureLoad.h"
#include "Debug.h"
#include "ElementGroup.h"
#include "PropertyTypeEnum.h"
#include "TimeFunction.h"


extern List<TimeFunction> timeFunction;


using namespace std;



Element2D3nodedPressureLoad::Element2D3nodedPressureLoad(void)
{
  if (debug) cout << " constructor Element2D3nodedPressureLoad\n\n";

  return;
}




Element2D3nodedPressureLoad::~Element2D3nodedPressureLoad()
{
  if (debug) cout << " destructor Element2D3nodedPressureLoad\n\n";

  return;
}







bool Element2D3nodedPressureLoad::forDomainType(int domType)
{
  switch (domType)
  {
    case MESH:      return true;
		
    case SOLID:     return true;
		
    default:        return false;
  }
}








int Element2D3nodedPressureLoad::finiteStrain(void)
{
  return roundToInt( ((ElementGroup*)elemGrp)->elemProp[ELEMENTTYPE]->data[2]);
}










int Element2D3nodedPressureLoad::calcStiffnessAndResidual(void)
{
  ElementGroup *eG = (ElementGroup*) elemGrp;
	
  double *s      = eG->dom->s,
         *p      = eG->dom->p,
         *xl     = eG->dom->xl,
         *ul     = eG->dom->ul,
         *elmDat = &(eG->elemProp[ELEMENTTYPE]->data[0]),
         press   = elmDat[1],
         x1[2], x2[2], x3[2], fact;

  int    i, j, l, 
         tmFctn  = roundToInt(elmDat[0]);

  bool   finite  = (roundToInt(elmDat[2])!=0),
         axsym   = (roundToInt(elmDat[3])!=0);


  // the pressure is applied at t_n!!!!!!!
  // !!!!!!! for dynamic problems some coding 
  // is required!!!!!!!

  press *= timeFunction[tmFctn].prop;
 
  for (i=0; i<3; i++)
  {
    for (j=0; j<2; j++)
    {
      ul[i+i+j]    =  u(i+1,j+1);
      xl[i+i+j]    =  x(i+1,j+1);
      xl[i+i+j+6]  = xn(i+1,j+1);
      xl[i+i+j+12] = x0(i+1,j+1);
    }
  }

  x1[0] = xl[12];
  x1[1] = xl[13];
  x2[0] = xl[14];
  x2[1] = xl[15];
  x3[0] = xl[16];
  x3[1] = xl[17];

  if (finite)
  {
    x1[0] += ul[0];
    x1[1] += ul[1];
    x2[0] += ul[2];
    x2[1] += ul[3];
    x3[0] += ul[4];
    x3[1] += ul[5];
  }

  for (i=0; i<36; i++) s[i] = 0.;

  fact = - press / 6.;

  p[0] = + fact * (3.*x1[1] - 4.*x2[1] + x3[1]);
  p[1] = - fact * (3.*x1[0] - 4.*x2[0] + x3[0]);

  p[4] = - fact * (3.*x3[1] - 4.*x2[1] + x1[1]);
  p[5] = + fact * (3.*x3[0] - 4.*x2[0] + x1[0]);

  fact *= 4.;

  p[2] = + fact * (x1[1] - x3[1]);
  p[3] = - fact * (x1[0] - x3[0]);

  if (!finite) return 0;

  fact *= - .25;

  s[ 6] = fact * 3.;
  s[18] = - fact * 4.;
  s[30] = fact;

  s[ 1] = - s[ 6];
  s[13] = - s[18];
  s[25] = - s[30];

  s[10] = s[25];
  s[22] = s[13];
  s[34] = s[ 1];

  s[ 5] = s[30];
  s[17] = s[18];
  s[29] = s[ 6];

  fact *= 4.;

  s[ 8] = fact;
  s[32] = - fact;

  s[ 3] = s[32];
  s[27] = s[ 8];

  return 0; 
}

