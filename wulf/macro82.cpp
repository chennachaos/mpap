
#include "Macro.h"
#include "FunctionsProgram.h"
#include "Plot.h"


extern Plot plot;


using namespace std;


int macro82(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "naca";
    macro.type = "wulf";
    macro.what = "generate NACA aerofoil geometetry";

    macro.sensitivity[INTER]  = true;
    macro.sensitivity[BATCH]  = true;
    macro.sensitivity[NOPROJ] = true;
    macro.sensitivity[PRE]    = true;
   
    macro.db.addList(COLOURS_BLUE);
 
    macro.db.addTextField("         4 digits: ",2412);
    macro.db.addTextField("    chord length = ",1.);
    macro.db.addTextField(" angle of attack = ",0.);
    macro.db.addTextField("number of points = ",30);
 
    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int    col   = roundToInt(macro.p[0]),
         naca  = roundToInt(macro.p[1]),
         n     = roundToInt(macro.p[4]*.5),
         i0, i1, i, d;
  double c     = macro.p[2],
         alpha = macro.p[3],
         t, p, m,
         xi, xi2, xi3, xi4, dd,
         x, y, yc, dyc,
         sn, cs, 
         xmx[2], xmn[2],
         *X = new double [n*4], dpt;
  char   tmp[10], tmp2[10];

  sprintf(tmp,"%d",naca);
  //cout << "|" << tmp << "|\n\n";
  i0 = strlen(tmp);
  //cout << i0;
  d  = 4 - i0;
  for (i=4; i>=i0; i--) tmp[i] = tmp[i-d];
  for (i=0; i<d; i++) tmp[i] = '0';
  //cout << "|" << tmp << "|\n\n";

  tmp2[0] = tmp[0]; tmp2[1] = '\0';
  if (!scanInt(tmp2,&d,false)) goto fail; m = .01 * (double) d;

  tmp2[0] = tmp[1];
  if (!scanInt(tmp2,&d,false)) goto fail; p = .1  * (double) d;

  tmp2[0] = tmp[2]; tmp2[1] = tmp[3]; tmp2[2] = '\0';
  if (!scanInt(tmp2,&d,false)) goto fail; t = .01 * (double) d;


  COUT << "    NACA " << tmp << "\n";
  COUT << "max. camber             : " << m << " c\n";
  COUT << "location of max. camber : " << p << " c\n";
  COUT << "max thickness in        : " << t << " c\n\n";

  dd = 3.14159265358979323846 / (double) n;

  for (i=0; i<=n; i++)
  {
    xi  = .5 * (1. + cos(i*dd));
    xi2 = xi*xi;
    xi3 = xi*xi2;
    xi4 = xi*xi3;

    y = 5.*t*c*(.2969*sqrt(xi)-.1281*xi-.3516*xi2+.2843*xi3-.1015*xi4);
 
    if (xi < p) 
    {
      yc  = m*c/(p*p) * xi*(p+p-xi);
      dyc = m*c/(p*p) * (p+p-xi-xi);
    }
    else
    {
      yc  = m * c/((1.-p)*(1.-p)) * (1.-xi)*(1.+xi-p-p);
      dyc = m * c/((1.-p)*(1.-p)) * (p+p-xi-xi);
    }

    cs = 1./sqrt(1.+dyc*dyc); 
    sn = cs * dyc;
 
    X[i+i]   = xi*c - y*sn; // xu
    X[i+i+1] = yc   + y*cs; // yu

    if (i != 0)
    {
      X[4*n-i-i-0] = xi*c + y*sn; // xl
      X[4*n-i-i+1] = yc   - y*cs; // yl
    }
  }

  alpha *= -0.017453293;
  sn = sin(alpha);
  cs = cos(alpha);

  for (i=0; i<n+n; i++)
  {
    x = X[i+i];
    y = X[i+i+1];

    X[i+i]   = c*.25 + (x-c*.25)*cs - y*sn;
    X[i+i+1] =         (x-c*.25)*sn + y*cs;
  }

  if (!plot)
  {
    xmn[0] = 0.;  xmn[1] = X[1];
    xmx[0] = c;   xmx[1] = X[1];

    for (i=1; i<n+n; i++) 
    {
      xmn[1] = min(xmn[1],X[i+i+1]);
      xmx[1] = max(xmx[1],X[i+i+1]);
    }
    plot.fit(xmn,xmx,2);
  }
  plot.setColour(col-1);
  
  dpt = (plot.dAct[0] + plot.dAct[1]) * .0035;
 
  for (i0=0; i0<n+n; i0++)
  {
    printf("    %3d %15.8f %15.8f\n",i0+1,X[i0+i0],X[i0+i0+1]);

    i1 = i0 + 1; if (i1 == n+n) i1 = 0;

    plot.point(X+i0+i0,dpt);

    plot.line(X+i0+i0,X+i1+i1);
  }
  cout << "\n";

  delete [] X;

  return 0;

  fail:

  delete [] X;

  prgWarning(1,"macro82","invalid NACA code!");

  return 0;  

//--------------------------------------------------------------------------------------------------
}

