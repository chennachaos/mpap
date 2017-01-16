
#include "Element1D2nodedPipeFlowST.h"
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



Element1D2nodedPipeFlowST::Element1D2nodedPipeFlowST(void)
{
  if (debug) cout << " constructor Element1D2nodedPipeFlowST\n\n";

  return;
}




Element1D2nodedPipeFlowST::~Element1D2nodedPipeFlowST()
{
  if (debug) cout << " destructor Element1D2nodedPipeFlowST\n\n";

  return;
}







bool Element1D2nodedPipeFlowST::forDomainType(int domType)
{
  switch (domType)
  {
    case MESH:      return true;
		
    case FEM1D:     return true;
		
    default:        return false;
  }
}




double power(double, int);


int Element1D2nodedPipeFlowST::calcStiffnessAndResidual(void)
{
/*  ElementGroup *eG = (ElementGroup*) elemGrp;

  int i, l;
  
  double *s      = eG->dom->s,
         *p      = eG->dom->p,
         *elmDat = &(eG->elemProp[ELEMENTTYPE]->data[0]),
	 rho, mu, h0, k, dt, dx,
         xx, tt, xi, eta, tau, dv, f1, f2, e, ex, ey, et,
         N1p, N1px, N1pt,
         N2p, N2px, N2pt,
         N1n, N1nx, N1nt,
         N2n, N2nx, N2nt,
         U, Ux, Ut,
         h, hx, ht,
         a, ax, ay, at, ah,
         uu, ux, uy, ut,
         vv, vy,
         px,
         U1p = un(1,1),
         U2p = un(2,1),
         U1n =  u(1,1),
         U2n =  u(2,1),
         h1p = un(1,2),
         h2p = un(2,2),
         h1n =  u(1,2),
         h2n =  u(2,2);


  rho  = elmDat[0];
  mu   = elmDat[1];
  R0   = elmDat[2];
  k    = elmDat[3];

  dt   = mpapTime.dt;
  dx   = x0(2,1) - x0(1,1);
 
  if (dx < 1.e-12) 
    prgError(1,"Element1D2nodedPipeFlowST::calcStiffnessAndResidual","dx < 0!");
 
  xx = 1. / dx;
  tt = 1. / dt;

  for (i=0; i<6;  i++) p[i] = 0.;
  for (i=0; i<36; i++) s[i] = 0.;

  for (l=0; l<8; l++)
  {
    // get Gauss point coordinates and weight

    f1 = .25;
    f2 = .75;

    switch (l)
    {
      case 0: xi = f1; eta = f1; tau = f1; break;
      case 1: xi = f2; eta = f1; tau = f1; break;
      case 2: xi = f2; eta = f2; tau = f1; break;
      case 3: xi = f1; eta = f2; tau = f1; break;
      case 4: xi = f1; eta = f1; tau = f2; break;  
      case 5: xi = f2; eta = f1; tau = f2; break;
      case 6: xi = f2; eta = f2; tau = f2; break;  
      case 7: xi = f1; eta = f2; tau = f2; break;
    }

    // shape functions

    Np   = 1. - tau;
    Npt  = - 1. / dt;

    Nn   = tau;
    Nnt  = 1. / dt;

    N1p  =   (1. - xi) * Np;
    N1px =             - Np / dx;
    N2pt = - (1. - xi) * Npt;

    N2p  =   xi * Np;
    N2px =        Np / dx;
    N2pt = - xi * Npt;

    N1n  =   (1. - xi) * Nn;
    N1nx =             - Nn / dx;
    N2nt =   (1. - xi) * Nnt;

    N2n  = xi * Nn;
    N2nx =      Nn / dx;
    N2nt = xi * Nnt;

    // radius R, axial and radial velocities U and V, pressure P and derivatives

    R    = N1p  * un(1,1) + N2p  * un(2,1) + N1n  * u(1,1) + N2n  * u(2,1) + R0;
    Rx   = N1px * un(1,1) + N2px * un(2,1) + N1nx * u(1,1) + N2nx * u(2,1);
    Rt   = N1pt * un(1,1) + N2pt * un(2,1) + N1nt * u(1,1) + N2nt * u(2,1);

    r    = R * eta;

    fct  = 1. - eta * eta;

    R1 

    fctx = (- eta - eta) * 

    U    = (N1p  * un(1,2) + N2p  * un(2,2) + N1n  * u(1,2) + N2n  * u(2,2)) * fct;
    Ux   = (N1px * un(1,2) + N2px * un(2,2) + N1nx * u(1,2) + N2nx * u(2,2)) * fct + U / fct * fctx;
    Ut   = (N1pt * un(1,2) + N2pt * un(2,2) + N1nt * u(1,2) + N2nt * u(2,2)) * fct;
    Ur   = U / fct * ;

    V    = Rt * eta;
    Vr   = ;
                                                                        
    P    = N1p  * un(1,3) + N2p  * un(2,3) + N1n  * u(1,3) + N2n  * u(2,3);
    Px  = N1px * un(1,3) + N2px * un(2,3) + N1nx * u(1,3) + N2nx * u(2,3);

    // standard Galerkin: Navier-Stokes in axial direction
    
    r[0] =  ( N1n * (rho * (Ut + V * Ur + U * Ux) + Px)
              + mu * (dUr1 * Ur + dUx1 * Ux - N1n * Ur / r) ) * r



    // standard Galerkin: continuity equation


    // standard Galerkin: elastic pipe



  }
*/
  return 0;
}

