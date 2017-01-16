
#include "Element1D2nodedFlexibleWing.h"
#include "Debug.h"
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



Element1D2nodedFlexibleWing::Element1D2nodedFlexibleWing(void)
{
  if (debug) cout << " constructor Element1D2nodedFlexibleWing\n\n";

  return;
}




Element1D2nodedFlexibleWing::~Element1D2nodedFlexibleWing()
{
  if (debug) cout << " destructor Element1D2nodedFlexibleWing\n\n";

  return;
}







bool Element1D2nodedFlexibleWing::forDomainType(int domType)
{
  switch (domType)
  {
    case MESH:      return true;
		
    case FEM1D:     return true;
		
    default:        return false;
  }
}







int Element1D2nodedFlexibleWing::calcStiffnessAndResidual(void)
{
  ElementGroup *eG = (ElementGroup*) elemGrp;
	
  double *s      = eG->dom->s,
         *p      = eG->dom->p,
         *elmDat = &(eG->elemProp[ELEMENTTYPE]->data[0]),
	 a, e, c, C, S, l, GJ, EI, a0, q;

  a   = elmDat[0];
  c   = elmDat[1];
  e   = elmDat[2];
  C   = cos(elmDat[3]/180.*3.1415927);
  S   = sin(elmDat[3]/180.*3.1415927);
  GJ  = elmDat[4];
  EI  = elmDat[5];
  q   = elmDat[6];
  a0  = elmDat[7];

  l   = x0(2,1) - x0(1,1);
 
  if (l < 1.e-12) 
    prgError(1,"Element1D2nodedFlexibleWing::calcStiffnessAndResidual","l < 0!");

  s[ 0] = 30*(24*EI + a*c*power(l, 3)*q*S);
  s[ 1] = 6*l*(60*EI + a*c*power(l, 3)*q*S);
  s[ 2] = 30*a*power(c, 2)*e*power(l, 3)*q*S;
  s[ 3] = 30*(-24*EI + a*c*power(l, 3)*q*S);
  s[ 4] = 360*EI*l - 6*a*c*power(l, 4)*q*S;
  s[ 5] = 30*a*power(c, 2)*e*power(l, 3)*q*S;
          
  s[ 6] = 360*EI*l - 6*a*c*power(l, 4)*q*S;
  s[ 7] = 240*EI*power(l, 2);
  s[ 8] = -5*a*power(c, 2)*e*power(l, 4)*q*S;
  s[ 9] = 6*l*(-60*EI + a*c*power(l, 3)*q*S);
  s[10] = 120*EI*power(l, 2) - a*c*power(l, 5)*q*S;
  s[11] = 5*a*power(c, 2)*e*power(l, 4)*q*S;
          
  s[12] = 21*a*c*C*power(l, 4)*q;
  s[13] = 3*a*c*C*power(l, 5)*q;
  s[14] = 20*power(l, 2)*(-3*GJ + a*power(c, 2)*C*e*power(l, 2)*q);
  s[15] = 9*a*c*C*power(l, 4)*q;
  s[16] = -2*a*c*C*power(l, 5)*q;
  s[17] = 10*power(l, 2)*(6*GJ + a*power(c, 2)*C*e*power(l, 2)*q);
          
  s[18] = -30*(24*EI + a*c*power(l, 3)*q*S);
  s[19] = -6*l*(60*EI + a*c*power(l, 3)*q*S);
  s[20] = -30*a*power(c, 2)*e*power(l, 3)*q*S;
  s[21] = 720*EI - 30*a*c*power(l, 3)*q*S;
  s[22] = 6*l*(-60*EI + a*c*power(l, 3)*q*S);
  s[23] = -30*a*power(c, 2)*e*power(l, 3)*q*S;
          
  s[24] = 6*l*(60*EI + a*c*power(l, 3)*q*S);
  s[25] = 120*EI*power(l, 2) + a*c*power(l, 5)*q*S;
  s[26] = 5*a*power(c, 2)*e*power(l, 4)*q*S;
  s[27] = -6*l*(60*EI + a*c*power(l, 3)*q*S);
  s[28] = 240*EI*power(l, 2);
  s[29] = -5*a*power(c, 2)*e*power(l, 4)*q*S;
          
  s[30] = 9*a*c*C*power(l, 4)*q;
  s[31] = 2*a*c*C*power(l, 5)*q;
  s[32] = 10*power(l, 2)*(6*GJ + a*power(c, 2)*C*e*power(l, 2)*q);
  s[33] = 21*a*c*C*power(l, 4)*q;
  s[34] = -3*a*c*C*power(l, 5)*q;
  s[35] = 20*power(l, 2)*(-3*GJ + a*power(c, 2)*C*e*power(l, 2)*q);

  p[0]  = -a0*30*a*c*power(l, 4)*q;
  p[1]  = -a0*5*a*c*power(l, 5)*q;
  p[2]  = -a0*30*a*power(c, 2)*e*power(l, 4)*q;
  p[3]  = -a0*30*a*c*power(l, 4)*q;
  p[4]  = -a0*(-5*a*c*power(l, 5)*q);
  p[5]  = -a0*30*a*power(c, 2)*e*power(l, 4)*q;

  int i, j, k;

  //prgPrintSimpleMatrix(s,6,6,12,5,true,1,false);
  //prgPrintSimpleMatrix(p,1,6,12,5,true,1,false);

  for (i=0; i<6; i++)
    for (j=0; j<2; j++)
      for (k=0; k<3; k++)
        p[i] -= s[(j*3+k)*6+i] * u(j+1,k+1);

  //prgPrintSimpleMatrix(p,1,6,12,5,true,1,false);

  return 0;
}

