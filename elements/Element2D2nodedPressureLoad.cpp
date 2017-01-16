
#include "Element2D2nodedPressureLoad.h"
#include "Debug.h"
#include "ElementGroup.h"
#include "PropertyTypeEnum.h"
#include "TimeFunction.h"


extern List<TimeFunction> timeFunction;


using namespace std;



Element2D2nodedPressureLoad::Element2D2nodedPressureLoad(void)
{
  if (debug) cout << " constructor Element2D2nodedPressureLoad\n\n";

  return;
}




Element2D2nodedPressureLoad::~Element2D2nodedPressureLoad()
{
  if (debug) cout << " destructor Element2D2nodedPressureLoad\n\n";

  return;
}







bool Element2D2nodedPressureLoad::forDomainType(int domType)
{
  switch (domType)
  {
    case MESH:      return true;
		
    case SOLID:     return true;
		
    default:        return false;
  }
}








int Element2D2nodedPressureLoad::finiteStrain(void)
{
  return roundToInt( ((ElementGroup*)elemGrp)->elemProp[ELEMENTTYPE]->data[2]);
}










int Element2D2nodedPressureLoad::calcStiffnessAndResidual(void)
{
  ElementGroup *eG = (ElementGroup*) elemGrp;
	
  double *s      = eG->dom->s,
         *p      = eG->dom->p,
         *xl     = eG->dom->xl,
         *ul     = eG->dom->ul,
         *elmDat = &(eG->elemProp[ELEMENTTYPE]->data[0]),
         press   = elmDat[1],
         dx, dy, sn, cs, fr;
 
  int    i, j,
         tmFctn  = roundToInt(elmDat[0]);

  bool   finite  = (roundToInt(elmDat[2])!=0),
         axsym   = (roundToInt(elmDat[3])!=0);


  // the pressure is applied at t_n!!!!!!!
  // !!!!!!! for dynamic problems some coding 
  // is required!!!!!!!

  press *= timeFunction[tmFctn].prop;
 
  for (i=0; i<2; i++)
  {
    for (j=0; j<2; j++)
    {
      ul[i+i+j]   =  u(i+1,j+1);
      xl[i+i+j]   =  x(i+1,j+1);
      xl[i+i+j+4] = xn(i+1,j+1);
      xl[i+i+j+8] = x0(i+1,j+1);
    }
  }

  dx = xl[8+2] - xl[8+0];
  dy = xl[8+3] - xl[8+1];

  if (finite) { dx += ul[2] - ul[0]; dy += ul[3] - ul[1]; }

  sn = dy;
  cs = dx;
  fr = press * .5;

  p[ 0] =   fr * sn; 
  p[ 1] = - fr * cs;
  p[ 2] = p[0];
  p[ 3] = p[1];

  if (!finite) { for (i=0; i<16; i++) s[i] = 0.; return 0; }

  s[ 0] = 0.;
  s[ 4] = fr;
  s[ 8] = 0.;
  s[12] = -fr;
  s[ 1] = -fr;
  s[ 5] = 0.;
  s[ 9] = fr;
  s[13] = 0.;
  s[ 2] = s[ 0];
  s[ 6] = s[ 4];
  s[10] = s[ 8];
  s[14] = s[12];
  s[ 3] = s[ 1];
  s[ 7] = s[ 5];
  s[11] = s[ 9];
  s[15] = s[13];

  return 0;
}

