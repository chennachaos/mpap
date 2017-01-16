
#include "Element2D2nodedFreeSurface.h"
#include "Plot.h"
#include "ElementGroup.h"
#include "FunctionsElement.h"
#include "FunctionsProgram.h"
#include "PropertyTypeEnum.h"
#include "FreeSurface.h"
#include "TimeFunction.h"


extern Plot plot;
extern List<TimeFunction> timeFunction;


using namespace std;



Element2D2nodedFreeSurface::Element2D2nodedFreeSurface(void)
{
  return;
}




Element2D2nodedFreeSurface::~Element2D2nodedFreeSurface()
{
  return;
}







bool Element2D2nodedFreeSurface::forDomainType(int domType)
{
  switch (domType)
  {
    case FREESURFACE: return true;
		
    default:          return false;
  }
}







int Element2D2nodedFreeSurface::nGaussPoints(void)
{
  if (((ElementGroup*)elemGrp)->elemProp[ELEMENTTYPE]->data.n < 1)

    prgError(1,"Element2D2nodedFreeSurface::nGaussPoints","element type data empty!");
	
  int ngp = roundToInt(((ElementGroup*)elemGrp)->elemProp[ELEMENTTYPE]->data[0]);

  if (ngp < 1 || ngp > 3) 
	  
    prgError(2,"Element2D2nodedFreeSurface::nGaussPoints","invalid number of Gauss points!");

  return ngp;
}






int Element2D2nodedFreeSurface::calcStiffnessAndResidual(void)
{
  ElementGroup *eG = (ElementGroup*) elemGrp;

  double *s = eG->dom->s,
         *p = eG->dom->p;

  int i;

  // maintain relative nodal spacing on free surface
  //   (elastic rope analogy)

  double dx0[2], dx[2], l2, l02, eps;

  dx0[0] = x0(2,1) - x0(1,1);
  dx0[1] = x0(2,2) - x0(1,2);

  dx[0] = x(2,1) - x(1,1);
  dx[1] = x(2,2) - x(1,2);

  l02 = dx0[0]*dx0[0] + dx0[1]*dx0[1];

  l2  = dx[0]*dx[0] + dx[1]*dx[1];

  eps = sqrt(l2/l02) - 1.;

  p[0] = eps;
  p[2] = -eps;

  s[ 0] = 0.;
  s[ 4] = 0.;
  s[ 8] = 0.;
  s[12] = 0.;
  s[16] = dx[0] / ((eps + 1.) * l02);
  s[20] = dx[1] / ((eps + 1.) * l02);
  s[24] = - s[16];
  s[28] = - s[20];

  s[ 2] = 0.;
  s[ 6] = 0.;
  s[10] = 0.;
  s[14] = 0.;
  s[18] = - s[16];
  s[22] = - s[20];
  s[26] = - s[18];
  s[30] = - s[22];

  // no fluid flow across free surface

  double ul[4], vl[4], xl[4], n[2], m[2], 
         td6  = eG->dom->td[5],
         td7  = eG->dom->td[6],
         td21 = eG->dom->td[20],
         fact1 = td7 * td21;

  ul[0] = td7 * u(1,1) + (1.-td7) * un(1,1);
  ul[1] = td7 * u(1,2) + (1.-td7) * un(1,2);
  ul[2] = td7 * u(2,1) + (1.-td7) * un(2,1);
  ul[3] = td7 * u(2,2) + (1.-td7) * un(2,2);

  vl[0] = td7 * v(1,1) + (1.-td7) * vn(1,1);
  vl[1] = td7 * v(1,2) + (1.-td7) * vn(1,2);
  vl[2] = td7 * v(2,1) + (1.-td7) * vn(2,1);
  vl[3] = td7 * v(2,2) + (1.-td7) * vn(2,2);

  xl[0] = td6 * x(1,1) + (1.-td6) * xn(1,1);
  xl[1] = td6 * x(1,2) + (1.-td6) * xn(1,2);
  xl[2] = td6 * x(2,1) + (1.-td6) * xn(2,1);
  xl[3] = td6 * x(2,2) + (1.-td6) * xn(2,2);

  dx[0] = xl[2] - xl[0];
  dx[1] = xl[3] - xl[1];

  n[0] = - dx[1];
  n[1] = + dx[0];

  m[0] = (ul[0] - vl[0]) * n[0] + (ul[1] - vl[1]) * n[1];
  m[1] = (ul[2] - vl[2]) * n[0] + (ul[3] - vl[3]) * n[1];

  p[1] = - m[0] - m[0] - m[1];
  p[3] = - m[0] - m[1] - m[1];

  s[ 1] = (n[0] + n[0]) * td7;
  s[ 5] = (n[1] + n[1]) * td7;
  s[ 9] = n[0] * td7;
  s[13] = n[1] * td7;
  s[17] = - (n[0] + n[0]) * fact1 - (ul[1]+ul[1]-vl[1]-vl[1]+ul[3]-vl[3]) * td6;
  s[21] = - (n[1] + n[1]) * fact1 + (ul[0]+ul[0]-vl[0]-vl[0]+ul[2]-vl[2]) * td6;
  s[25] = - n[0] * fact1          + (ul[1]+ul[1]-vl[1]-vl[1]+ul[3]-vl[3]) * td6;
  s[29] = - n[1] * fact1          - (ul[0]+ul[0]-vl[0]-vl[0]+ul[2]-vl[2]) * td6;

  s[ 3] = n[0] * td7;
  s[ 7] = n[1] * td7;
  s[11] = (n[0] + n[0]) * td7;
  s[15] = (n[1] + n[1]) * td7;
  s[19] = - n[0] * fact1          - (ul[1]-vl[1]+ul[3]+ul[3]-vl[3]-vl[3]) * td6;;
  s[23] = - n[1] * fact1          + (ul[0]-vl[0]+ul[2]+ul[2]-vl[2]-vl[2]) * td6;;
  s[27] = - (n[0] + n[0]) * fact1 + (ul[1]-vl[1]+ul[3]+ul[3]-vl[3]-vl[3]) * td6;;
  s[31] = - (n[1] + n[1]) * fact1 - (ul[0]-vl[0]+ul[2]+ul[2]-vl[2]-vl[2]) * td6;;

  return 0;
}

