
#include "Element1D2nodedAdvectionDiffusion.h"
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



Element1D2nodedAdvectionDiffusion::Element1D2nodedAdvectionDiffusion(void)
{
  if (debug) cout << " constructor Element1D2nodedAdvectionDiffusion\n\n";

  return;
}




Element1D2nodedAdvectionDiffusion::~Element1D2nodedAdvectionDiffusion()
{
  if (debug) cout << " destructor Element1D2nodedAdvectionDiffusion\n\n";

  return;
}







bool Element1D2nodedAdvectionDiffusion::forDomainType(int domType)
{
  switch (domType)
  {
    case MESH:      return true;
		
    case FEM1D:     return true;
		
    default:        return false;
  }
}







int Element1D2nodedAdvectionDiffusion::calcStiffnessAndResidual(void)
{
  ElementGroup *eG = (ElementGroup*) elemGrp;
	
  double *s      = eG->dom->s,
         *p      = eG->dom->p,
         *elmDat = &(eG->elemProp[ELEMENTTYPE]->data[0]),
	 A1, A2, A3, A4,
	 a, gam, mu, dt, dx, tau;

/* 
   1D advection diffusion 
   ----------------------

   generalised midpoint rule

   see PhD p.95

*/
  
  a   = elmDat[0];
  mu  = elmDat[1];
  gam = elmDat[2];
  dt  = mpapTime.dt;
  dx  = x0(2,1) - x0(1,1);
  tau = 0.;  
 
  if (dx < 1.e-12) 
    prgError(1,"Element1D2nodedAdvectionDiffusion::calcStiffnessAndResidual","dx < 0!");
 
  A1 = 1./6.;
  A2 = tau * a * 0.5 / dx;
  A3 = gam * a * dt * 0.5 / dx;
  A4 = gam * dt * (mu + tau * a*a) / (dx*dx);
	  
  s[0] = A1 + A1 - A2 - A3 + A4;
  s[1] =      A1 + A2 - A3 - A4;
  s[2] =      A1 - A2 + A3 - A4;
  s[3] = A1 + A1 + A2 + A3 + A4;

  p[0] = - s[0] * u(1,1) - s[2] * u(2,1);
  p[1] = - s[1] * u(1,1) - s[3] * u(2,1);
  
  A1 = - A1;
  A2 = - A2;
  A3 = A3 / gam * (1. - gam);
  A4 = A4 / gam * (1. - gam);
  
  p[0] += - (A1 + A1 - A2 - A3 + A4) * un(1,1) - (A1 - A2 + A3 - A4) * un(2,1);
  p[1] += - (A1 + A2 - A3 - A4) * un(1,1) - (A1 + A1 + A2 + A3 + A4) * un(2,1);

  return 0;
}

