
#include "Element1D2nodedSchroedingerCN.h"
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



Element1D2nodedSchroedingerCN::Element1D2nodedSchroedingerCN(void)
{
  if (debug) cout << " constructor Element1D2nodedSchroedingerCN\n\n";

  return;
}




Element1D2nodedSchroedingerCN::~Element1D2nodedSchroedingerCN()
{
  if (debug) cout << " destructor Element1D2nodedSchroedingerCN\n\n";

  return;
}







bool Element1D2nodedSchroedingerCN::forDomainType(int domType)
{
  switch (domType)
  {
    case MESH:      return true;
		
    case FEM1D:     return true;
		
    default:        return false;
  }
}




double power(double, int);


int Element1D2nodedSchroedingerCN::calcStiffnessAndResidual(void)
{
  ElementGroup *eG = (ElementGroup*) elemGrp;

  int i;
  
  double *s      = eG->dom->s,
         *p      = eG->dom->p,
         *elmDat = &(eG->elemProp[ELEMENTTYPE]->data[0]),
	 alp, bet,
	 h, m, fact, dt, dx;

/* 
   1D SchroedingerCN equation
   --------------------------

   Crank-Nicolson (trapezoidal rule)

   Avtar's MRes

*/
  
  h    = elmDat[0];
  m    = elmDat[1];
  fact = elmDat[2];
  dt   = mpapTime.dt;
  dx   = x0(2,1) - x0(1,1);
 
  if (dx < 1.e-12) 
    prgError(1,"Element1D2nodedSchroedingerCN::calcStiffnessAndResidual","dx < 0!");
 
  alp = dt * h / (4. * m * dx);
  bet = dx / 6.;

  s[ 0] = bet+bet;
  s[ 1] = alp;
  s[ 2] = bet;
  s[ 3] = - alp;
  
  s[ 4] = - alp;
  s[ 5] = bet+bet;
  s[ 6] = alp;
  s[ 7] = bet;
  
  s[ 8] = bet;
  s[ 9] = - alp;
  s[10] = bet+bet;
  s[11] = alp;
  
  s[12] = alp;
  s[13] = bet;
  s[14] = - alp;
  s[15] = bet+bet;
 
  p[0] = - s[0] * u(1,1) - s[4] * u(1,2) - s[ 8] * u(2,1) - s[12] * u(2,2);
  p[1] = - s[1] * u(1,1) - s[5] * u(1,2) - s[ 9] * u(2,1) - s[13] * u(2,2);
  p[2] = - s[2] * u(1,1) - s[6] * u(1,2) - s[10] * u(2,1) - s[14] * u(2,2);
  p[3] = - s[3] * u(1,1) - s[7] * u(1,2) - s[11] * u(2,1) - s[15] * u(2,2);
  
  p[0] += + s[0] * un(1,1) - s[4] * un(1,2) + s[ 8] * un(2,1) - s[12] * un(2,2);
  p[1] += - s[1] * un(1,1) + s[5] * un(1,2) - s[ 9] * un(2,1) + s[13] * un(2,2);
  p[2] += + s[2] * un(1,1) - s[6] * un(1,2) + s[10] * un(2,1) - s[14] * un(2,2);
  p[3] += - s[3] * un(1,1) + s[7] * un(1,2) - s[11] * un(2,1) + s[15] * un(2,2);


  double a1  =  u(1,1),
	 b1  =  u(1,2),
	 a2  =  u(2,1),
         b2  =  u(2,2),
	 a1n = un(1,1),
	 b1n = un(1,2),
	 a2n = un(2,1),
         b2n = un(2,2);
 

  p[0] += -((2*power(a2,2)*b1 + 4*a2*a2n*b1 + 2*power(a2n,2)*b1 + 12*power(b1,3) + 
			          2*power(a2,2)*b1n + 4*a2*a2n*b1n + 2*power(a2n,2)*b1n + 36*power(b1,2)*b1n + 
				          36*b1*power(b1n,2) + 12*power(b1n,3) + 3*power(a2,2)*b2 + 6*a2*a2n*b2 + 
					          3*power(a2n,2)*b2 + 9*power(b1,2)*b2 + 18*b1*b1n*b2 + 9*power(b1n,2)*b2 + 
						          6*b1*power(b2,2) + 6*b1n*power(b2,2) + 3*power(b2,3) + 
							          3*(power(a2 + a2n,2) + 3*power(b1 + b1n,2) + 4*(b1 + b1n)*b2 + 3*power(b2,2))*
								           b2n + 3*(2*(b1 + b1n) + 3*b2)*power(b2n,2) + 3*power(b2n,3) + 
									           3*power(a1,2)*(4*b1 + 4*b1n + b2 + b2n) + 
										           3*power(a1n,2)*(4*b1 + 4*b1n + b2 + b2n) + 
											           2*a1n*(a2 + a2n)*(3*b1 + 3*b1n + 2*(b2 + b2n)) + 
												           a1*(6*a1n*(4*b1 + 4*b1n + b2 + b2n) + 2*(a2 + a2n)*(3*b1 + 3*b1n + 2*(b2 + b2n)))
													           )*dt*dx*fact)/480.;
	  
  p[1] += ((12*power(a1,3) + 12*power(a1n,3) + 9*power(a1n,2)*(a2 + a2n) + 
			         9*power(a1,2)*(4*a1n + a2 + a2n) + 
				        2*a1n*(3*power(a2 + a2n,2) + 6*power(b1,2) + 6*power(b1n,2) + 3*b1n*(b2 + b2n) + 
						          power(b2 + b2n,2) + 3*b1*(4*b1n + b2 + b2n)) + 
					       2*a1*(3*(6*power(a1n,2) + 3*a1n*(a2 + a2n) + power(a2 + a2n,2)) + 6*power(b1,2) + 
						                 6*power(b1n,2) + 3*b1n*(b2 + b2n) + power(b2 + b2n,2) + 3*b1*(4*b1n + b2 + b2n)
								           ) + (a2 + a2n)*(3*power(a2 + a2n,2) + 3*power(b1,2) + 3*power(b1n,2) + 
										             4*b1n*(b2 + b2n) + 3*power(b2 + b2n,2) + b1*(6*b1n + 4*(b2 + b2n))))*dt*dx*fact)/
	     480.;
  
  p[2] += -((power(a1,2)*(3*b1 + 3*b1n + 2*(b2 + b2n)) + 
			          power(a1n,2)*(3*b1 + 3*b1n + 2*(b2 + b2n)) + 
				          2*a1n*(a2 + a2n)*(2*b1 + 2*b1n + 3*(b2 + b2n)) + 
					          3*(power(b1,3) + 3*power(b1,2)*b1n + 3*b1*power(b1n,2) + power(b1n,3) + 
							             2*power(b1,2)*b2 + 4*b1*b1n*b2 + 2*power(b1n,2)*b2 + 3*b1*power(b2,2) + 
								                3*b1n*power(b2,2) + 4*power(b2,3) + 
										           2*(power(b1 + b1n,2) + 3*(b1 + b1n)*b2 + 6*power(b2,2))*b2n + 
											              3*(b1 + b1n + 4*b2)*power(b2n,2) + 4*power(b2n,3) + 
												                 power(a2,2)*(b1 + b1n + 4*(b2 + b2n)) + 2*a2*a2n*(b1 + b1n + 4*(b2 + b2n)) + 
														            power(a2n,2)*(b1 + b1n + 4*(b2 + b2n))) + 
						          a1*(2*(a2 + a2n)*(2*b1 + 2*b1n + 3*(b2 + b2n)) + 
								             a1n*(6*b1 + 6*b1n + 4*(b2 + b2n))))*dt*dx*fact)/480.;
  
  p[3] += ((3*power(a1,3) + 3*power(a1n,3) + 6*power(a1n,2)*(a2 + a2n) + 
			         3*power(a1,2)*(3*a1n + 2*(a2 + a2n)) + 
				        2*(a2 + a2n)*(6*power(a2 + a2n,2) + power(b1,2) + power(b1n,2) + 
						          3*b1n*(b2 + b2n) + 6*power(b2 + b2n,2) + b1*(2*b1n + 3*(b2 + b2n))) + 
					       a1n*(9*power(a2 + a2n,2) + 3*power(b1,2) + 3*power(b1n,2) + 4*b1n*(b2 + b2n) + 
						                 3*power(b2 + b2n,2) + b1*(6*b1n + 4*(b2 + b2n))) + 
					              a1*(3*(3*power(a1n,2) + 4*a1n*(a2 + a2n) + 3*power(a2 + a2n,2)) + 3*power(b1,2) + 
							                3*power(b1n,2) + 4*b1n*(b2 + b2n) + 3*power(b2 + b2n,2) + 
									          b1*(6*b1n + 4*(b2 + b2n))))*dt*dx*fact)/480.;


  s[0] -= -((3*a1*(4*b1 + 4*b1n + b2 + b2n) + 3*a1n*(4*b1 + 4*b1n + b2 + b2n) + 
			          (a2 + a2n)*(3*b1 + 3*b1n + 2*(b2 + b2n)))*dt*dx*fact)/240.;
  
  s[1] -= ((3*(6*power(a1,2) + 6*power(a1n,2) + 3*a1n*(a2 + a2n) + power(a2 + a2n,2) + 
				            3*a1*(4*a1n + a2 + a2n)) + 6*power(b1,2) + 6*power(b1n,2) + 3*b1n*(b2 + b2n) + 
			         power(b2 + b2n,2) + 3*b1*(4*b1n + b2 + b2n))*dt*dx*fact)/240.;
  
  s[2] -= -((a1*(3*b1 + 3*b1n + 2*(b2 + b2n)) + a1n*(3*b1 + 3*b1n + 2*(b2 + b2n)) + 
			          (a2 + a2n)*(2*b1 + 2*b1n + 3*(b2 + b2n)))*dt*dx*fact)/240.;
  
  s[3] -= ((3*(3*power(a1,2) + 3*power(a1n,2) + 4*a1n*(a2 + a2n) + 3*power(a2 + a2n,2) + 
				            a1*(6*a1n + 4*(a2 + a2n))) + 3*power(b1,2) + 3*power(b1n,2) + 
			         4*b1n*(b2 + b2n) + 3*power(b2 + b2n,2) + b1*(6*b1n + 4*(b2 + b2n)))*dt*dx*fact)/480.;

  
  s[4] -= -((6*power(a1,2) + 6*power(a1n,2) + 3*a1n*(a2 + a2n) + power(a2 + a2n,2) + 
			          3*a1*(4*a1n + a2 + a2n) + 3*
				           (6*power(b1,2) + 6*power(b1n,2) + 3*b1n*(b2 + b2n) + power(b2 + b2n,2) + 
					               3*b1*(4*b1n + b2 + b2n)))*dt*dx*fact)/240.;
  
  s[5] -= ((3*a1*(4*b1 + 4*b1n + b2 + b2n) + 3*a1n*(4*b1 + 4*b1n + b2 + b2n) + 
			         (a2 + a2n)*(3*b1 + 3*b1n + 2*(b2 + b2n)))*dt*dx*fact)/240.;
  
  s[6] -= -((3*power(a1,2) + 3*power(a1n,2) + 4*a1n*(a2 + a2n) + a1*(6*a1n + 4*(a2 + a2n)) + 
			          3*(power(a2 + a2n,2) + 3*power(b1,2) + 3*power(b1n,2) + 4*b1n*(b2 + b2n) + 
					             3*power(b2 + b2n,2) + b1*(6*b1n + 4*(b2 + b2n))))*dt*dx*fact)/480.;

  s[7] -= ((a1*(3*b1 + 3*b1n + 2*(b2 + b2n)) + a1n*(3*b1 + 3*b1n + 2*(b2 + b2n)) + 
			         (a2 + a2n)*(2*b1 + 2*b1n + 3*(b2 + b2n)))*dt*dx*fact)/240.;
  
  
  s[8] -= -((a1*(3*b1 + 3*b1n + 2*(b2 + b2n)) + a1n*(3*b1 + 3*b1n + 2*(b2 + b2n)) + 
			          (a2 + a2n)*(2*b1 + 2*b1n + 3*(b2 + b2n)))*dt*dx*fact)/240.;
  
  s[9] -= ((3*(3*power(a1,2) + 3*power(a1n,2) + 4*a1n*(a2 + a2n) + 3*power(a2 + a2n,2) + 
				            a1*(6*a1n + 4*(a2 + a2n))) + 3*power(b1,2) + 3*power(b1n,2) + 
			         4*b1n*(b2 + b2n) + 3*power(b2 + b2n,2) + b1*(6*b1n + 4*(b2 + b2n)))*dt*dx*fact)/480.;
  
  s[10]-= -((a1*(2*b1 + 2*b1n + 3*(b2 + b2n)) + a1n*(2*b1 + 2*b1n + 3*(b2 + b2n)) + 
			          3*(a2 + a2n)*(b1 + b1n + 4*(b2 + b2n)))*dt*dx*fact)/240.;
  
  s[11]-= ((3*(power(a1,2) + power(a1n,2) + 3*a1n*(a2 + a2n) + 6*power(a2 + a2n,2) + 
				            a1*(2*a1n + 3*(a2 + a2n))) + power(b1,2) + power(b1n,2) + 3*b1n*(b2 + b2n) + 
			         6*power(b2 + b2n,2) + b1*(2*b1n + 3*(b2 + b2n)))*dt*dx*fact)/240.;

  
  s[12]-= -((3*power(a1,2) + 3*power(a1n,2) + 4*a1n*(a2 + a2n) + a1*(6*a1n + 4*(a2 + a2n)) + 
			          3*(power(a2 + a2n,2) + 3*power(b1,2) + 3*power(b1n,2) + 4*b1n*(b2 + b2n) + 
					             3*power(b2 + b2n,2) + b1*(6*b1n + 4*(b2 + b2n))))*dt*dx*fact)/480.;
  
  s[13]-= ((a1*(3*b1 + 3*b1n + 2*(b2 + b2n)) + a1n*(3*b1 + 3*b1n + 2*(b2 + b2n)) + 
			         (a2 + a2n)*(2*b1 + 2*b1n + 3*(b2 + b2n)))*dt*dx*fact)/240.;
  
  s[14]-= -((power(a1,2) + power(a1n,2) + 3*a1n*(a2 + a2n) + a1*(2*a1n + 3*(a2 + a2n)) + 
			          3*(2*power(a2 + a2n,2) + power(b1,2) + power(b1n,2) + 3*b1n*(b2 + b2n) + 
					             6*power(b2 + b2n,2) + b1*(2*b1n + 3*(b2 + b2n))))*dt*dx*fact)/240.;
  
  s[15]-= ((a1*(2*b1 + 2*b1n + 3*(b2 + b2n)) + a1n*(2*b1 + 2*b1n + 3*(b2 + b2n)) + 
			         3*(a2 + a2n)*(b1 + b1n + 4*(b2 + b2n)))*dt*dx*fact)/240.;

  
  return 0;
}

