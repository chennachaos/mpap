
#include <iostream>
#include <cmath>

#include "FunctionsProgram.h"
#include "DAFunctionPlot.h"
#include "ColoursEnum.h"


using namespace std;


#define N 100



void optimizeFlightPath(void)
{
  DAFunctionPlot plot;
/*
  double R, H, T, 
         dRdg, dRdT, dRdg2, dRdT2,
         dHdg, dHdT, dHdg2, dHdT2,
         v, g, dv, dg, x, h,
         vn, gn, vn1, gn1, 
         r1, r2, dr1dv, dr1dg, dr2dv, dr2dg, 
         fact, tol, dt, S, rho, CD, m, grav;

  int n;

  tol  = 1.e-10;
  dt   = 0.003;
  rho  = 1.;
  CD   = 100.;
  S    = 1.;
  m    = 1.;
  grav = -1.;

  vn1  = 1.;
  gn1  = .1;

  h = 0.;
  x = 0.;

  plot.addXY(x,h);

  for (n=1; n<=N; n++)
  {
    vn = vn1;

    gn = gn1;

    while (1)
    {
      v     = (vn1 + vn) * .5;
           
      g     = (gn1 + gn) * .5;
           
      dv    = (vn1 - vn) / dt;
           
      dg    = (gn1 - gn) / dt;
           
      r1    = dv + rho * S * CD / (2. * m) * v * v + grav * sin(g);
           
      r2    = v * dg - grav * cos(g);

      cout << sqrt(r1*r1+r2*r2) << "\n";

      if (sqrt(r1*r1+r2*r2) < tol) break;
      
      dr1dv = 1. / dt + rho * S * CD / (2. * m) * v;

      dr1dg = grav * cos(g) * .5;

      dr2dv = dg * .5;

      dr2dg = v / dt + grav * sin(g) * .5;

      fact  = 1. / (dr1dv * dr2dg - dr1dg * dr2dv);

      vn1  -= fact * (r1 * dr2dg - dr1dg * r2);

      gn1  -= fact * (dr1dv * r2 - r1 * dr2dv);
    }

    //prgUpdateDisplay(true);

    cout << "\n";

    x += v * dt * cos(g);

    h += v * dt * sin(g);

    dxdg  = dxdg + dt * (dvdg * cos(g) - v * sin(g) * dgdg);

    dxdt  = dxdt + dt * (dvdt * cos(g) - v * sin(g) * dgdt) + v * cos(g);

    dxdgg = 

    dxdgt = 

    dxdtg = 

    dxdtt = 

    dhdg  = dhdg + dt * (dvdg * sin(g) + v * cos(g) * dgdg);

    dhdt  = dhdt + dt * (dvdt * sin(g) + v * cos(g) * dgdt) + v * sin(g);

    if (h < 0.) break;

    plot.addXY(x,h);
  }

  plot.setup(10,10,500,300,BLUE,YELLOW,BLUE,5);

  plot.draw();

  prgUpdateDisplay();
*/
  return;
}





