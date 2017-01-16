
#include "GeomSpline.h"
#include "Plot.h"
#include "MathVector.h"
#include "Geometry.h"
#include "MathGeom.h"



extern Plot plot;




GeomSpline::GeomSpline(void)
{
  return;
}




GeomSpline::~GeomSpline()
{
  return;
}







double GeomSpline::xmin(void)
{
  double xmn = ((GeomPoint*)(suppPnt[0]))->x[0];

  for (int i=1; i<suppPnt.n; i++) 
    if (xmn > ((GeomPoint*)(suppPnt[i]))->x[0]) xmn = ((GeomPoint*)(suppPnt[i]))->x[0];

  return xmn;
}








double GeomSpline::area(double *xp)
{
  double ar = 0., *x0, *x1;

  // this is based on the geometry points only (-> inaccurate)
  // and should be used only to check the sign of the area

  x0 = ((GeomPoint*)(suppPnt[0]))->x;

  for (int i=1; i<suppPnt.n; i++)
  {
    x1  = ((GeomPoint*)(suppPnt[i]))->x;

    ar += triangleArea2D(xp,x0,x1);

    x0 = x1;
  }

  return ar;
}







void GeomSpline::discretise(Domain *dom, int nd)
{
  char fct[] = "GeomSpline::discretise";

  ListInfinite<VectorFixed<double,5> > splDiv;

  List<GeomPoint> &point = ((Geometry*)dom)->point;

  int c, i, j, n, npt, mult;

  if (suppPnt.n < 1) return;

  double x1[2], x2[2], d, *X, minh, a1, a2, m = 0., hfact, fact;

  // find minimum element size

  minh = ((GeomPoint*)(suppPnt[0]))->h;
  for (i=1; i<suppPnt.n; i++)
    if (minh > ((GeomPoint*)(suppPnt[i]))->h) minh = ((GeomPoint*)(suppPnt[i]))->h;

  //cout << minh << " = minh\n\n";

  if      (nd == 3) { hfact = 1.0; mult = 1; }
  else if (nd == 4) { hfact = 1.0; mult = 2; }
  else if (nd == 6) { hfact = 0.5; mult = 2; }
  else if (nd == 8) { hfact = 0.5; mult = 4; }
  else if (nd == 9) { hfact = 0.5; mult = 4; }
  else prgError(1,fct,"invalid value of nd!");

  minh *= hfact;

  // generate subdivision of the spline

  calcTangent(0,tang1);

  x1[0] = ((GeomPoint*)suppPnt[0])->x[0];
  x1[1] = ((GeomPoint*)suppPnt[0])->x[1];

  X = splDiv[0].x;
  X[0] = x1[0];
  X[1] = x1[1];
  X[2] = hfact * ((GeomPoint*)suppPnt[0])->h;
  X[3] = 0.;
  X[4] = 0.;

  for (i=1; i<suppPnt.n; i++)   
  {
    calcTangent(i,tang2);

    x2[0] = ((GeomPoint*)suppPnt[i])->x[0];
    x2[1] = ((GeomPoint*)suppPnt[i])->x[1];

    d = sqrt((x2[0]-x1[0])*(x2[0]-x1[0]) + (x2[1]-x1[1])*(x2[1]-x1[1]));

    if (d > minh) subdivide(i,splDiv,hfact,minh,x1,x2,0.,1.);

    else
    {
      X = splDiv[splDiv.n].x;
      X[0] = x2[0];
      X[1] = x2[1];
      X[2] = hfact * ((GeomPoint*)suppPnt[i])->h;
      X[3] = d;
      X[4] = (double)i;
    }

    x1[0] = x2[0];
    x1[1] = x2[1];

    tang1[0] = tang2[0];
    tang1[1] = tang2[1];
  }

  // generate new geometry points

  for (i=1; i<splDiv.n; i++) m += 2. * splDiv[i][3] / (splDiv[i-1][2] + splDiv[i][2]);

  //n = mult; while ((double)n < m-.1) n += mult;

  n = max(1, roundToInt(m / ((double)mult)) * mult);

  //cout << m << "->" << n << "\n";

  fact = m / ((double)n);

  for (i=0; i<splDiv.n; i++) splDiv[i][2] *= fact;

  //for (i=0; i<splDiv.n; i++) cout << i << ":" << splDiv[i][3] << "," << splDiv[i][2] << "\n";

  splDiv[0][3] = 0.;
  for (i=1; i<splDiv.n; i++) 
    splDiv[i][3] = splDiv[i-1][3] + splDiv[i][3] * 2. / (splDiv[i-1][2] + splDiv[i][2]);

  Vector<void*> suppPntTmp;

  npt = point.n;
  c = 0;
  i = 1;

  while (++c != n)
  {
    while (i < splDiv.n && ((double)c)-.0001 > splDiv[i][3]) i++;

    a1   = (double)c - splDiv[i-1][3];
    a2   = splDiv[i][3] - (double)c;
    fact = 1. / (a1 + a2);
    a1  *= fact;
    a2  *= fact;

    point.add(new GeomPoint);

    fact = a1 * splDiv[i][4] + a2 * splDiv[i-1][4];

    j = (int)(floor(fact));

    fact -= floor(fact);

    calcTangent(j,  tang1);
    calcTangent(j+1,tang2);

    calcPoint  (j+1,fact, point[npt].x);
  
    point[npt].h = (a1 * splDiv[i][2] + a2 * splDiv[i-1][2]) / hfact;

    point[npt].transferFact.append(1. - fact);
    point[npt].transferFact.append(fact);

    point[npt].transferNode.append(((GeomPoint*)suppPnt[j  ])->dat1 + 1);
    point[npt].transferNode.append(((GeomPoint*)suppPnt[j+1])->dat1 + 1);
 
    suppPntTmp.append((void*)&(point[npt++]));
  }

  ((GeomPoint*)suppPnt[0])->transferNode.free();
  ((GeomPoint*)suppPnt[0])->transferFact.free();

  ((GeomPoint*)suppPnt[0])->transferNode.append(((GeomPoint*)suppPnt[0])->dat1 + 1);
  ((GeomPoint*)suppPnt[0])->transferFact.append(1.);

  ((GeomPoint*)suppPnt.lastCoeff())->transferNode.free();
  ((GeomPoint*)suppPnt.lastCoeff())->transferFact.free();

  ((GeomPoint*)suppPnt.lastCoeff())->transferNode.append(((GeomPoint*)suppPnt.lastCoeff())->dat1+1);
  ((GeomPoint*)suppPnt.lastCoeff())->transferFact.append(1.);

  // delete old geometry points

  for (i=1; i<suppPnt.n-1; i++) point.del((GeomPoint*)suppPnt[i]);

  suppPntTmp.append(suppPnt[suppPnt.n-1]);

  suppPnt.trunc(1);

  for (i=0; i<suppPntTmp.n; i++) suppPnt.append(suppPntTmp[i]); 

  return;
}





void GeomSpline::subdivide(int iSect, ListInfinite<VectorFixed<double,5> > &splDiv, 
                           double hfact, double minh, double *x1, double *x2, 
                           double t1, double t2)
{
  double x3[2], *X, t = .5 * (t1 + t2), d;

  calcPoint(iSect,t,x3);

  d = sqrt((x3[0]-x1[0])*(x3[0]-x1[0]) + (x3[1]-x1[1])*(x3[1]-x1[1]));

  if (d > minh) subdivide(iSect,splDiv,hfact,minh,x1,x3,t1,t);
  
  else
  {
    X = splDiv[splDiv.n].x;
    X[0] = x3[0];
    X[1] = x3[1];
    X[2] = hfact * (((GeomPoint*)suppPnt[iSect-1])->h * (1.-t)
                  + ((GeomPoint*)suppPnt[iSect  ])->h *   t   );
    X[3] = d;
    X[4] = (double)iSect - 1. + t;
  }

  d = sqrt((x2[0]-x3[0])*(x2[0]-x3[0]) + (x2[1]-x3[1])*(x2[1]-x3[1]));

  if (d > minh) subdivide(iSect,splDiv,hfact,minh,x3,x2,t,t2);

  else
  {
    X = splDiv[splDiv.n].x;
    X[0] = x2[0];
    X[1] = x2[1];
    X[2] = hfact * (((GeomPoint*)suppPnt[iSect-1])->h * (1.-t2)
                  + ((GeomPoint*)suppPnt[iSect  ])->h *   t2   );
    X[3] = d;
    X[4] = (double)iSect - 1. + t2;
  }

  return;
}








void GeomSpline::checkElemSizes(void)
{
  int i;

  double d = 0.;

  for (i=1; i<suppPnt.n; i++)
    d += sqrt(dist2(((GeomPoint*)(suppPnt[i]))->x,((GeomPoint*)(suppPnt[i-1]))->x,2));

  for (i=0; i<suppPnt.n; i++)
    ((GeomPoint*)(suppPnt[i]))->h = min(((GeomPoint*)(suppPnt[i]))->h,d);

  return;
}








void GeomSpline::draw(int splineId)
{
  int i, j, n = 100;

  if (suppPnt.n < 1) return;

  double x1[2], x2[2], t, dt = 1./ ((double) n);

  calcTangent(0,tang1);

  x1[0] = ((GeomPoint*)suppPnt[0])->x[0];
  x1[1] = ((GeomPoint*)suppPnt[0])->x[1];

  for (i=1; i<suppPnt.n; i++)   
  {
    calcTangent(i,tang2);

    t = 0.;

    for (j=1; j<n; j++)
    {
      t += dt;

      calcPoint(i,t,x2);

      plot.line(x1,x2);

      x1[0] = x2[0];
      x1[1] = x2[1];
    }
    x2[0] = ((GeomPoint*)suppPnt[i])->x[0];
    x2[1] = ((GeomPoint*)suppPnt[i])->x[1];

    plot.line(x1,x2);

    x1[0] = x2[0];
    x1[1] = x2[1];

    tang1[0] = tang2[0];
    tang1[1] = tang2[1];
  }

  return;
}





void GeomSpline::calcPoint(int iSect, double t, double *x)
{
  double *x1 = ((GeomPoint*)suppPnt[iSect-1])->x, 
         *x2 = ((GeomPoint*)suppPnt[iSect  ])->x,
         t2  = t * t, 
         t3  = t2 * t,
         h2  = t2+t2+t2 - t3-t3,
         h1  = 1. - h2,
         h4  = t3 - t2,
         h3  = h4 - t2 + t;

  x[0] = h1 * x1[0] + h2 * x2[0] + h3 * tang1[0] + h4 * tang2[0];
  x[1] = h1 * x1[1] + h2 * x2[1] + h3 * tang1[1] + h4 * tang2[1];

  return;
}





void GeomSpline::calcTangent(int iPnt, double *tang)
{
  if (suppPnt.n == 2) 
  {
    tang[0] = ((GeomPoint*)suppPnt[1])->x[0] - ((GeomPoint*)suppPnt[0])->x[0];
    tang[1] = ((GeomPoint*)suppPnt[1])->x[1] - ((GeomPoint*)suppPnt[0])->x[1];
    return;
  }

  int i = iPnt;
  if (iPnt <= 0)           i = 1;
  if (iPnt >= suppPnt.n-1) i = suppPnt.n - 2;

  double *x1 = ((GeomPoint*)suppPnt[i-1])->x, 
         *x2 = ((GeomPoint*)suppPnt[i  ])->x,
         *x3 = ((GeomPoint*)suppPnt[i+1])->x;

  tang[0] = (x3[0] - x1[0])*.45;
  tang[1] = (x3[1] - x1[1])*.45;

  if (iPnt == i) return;

  double alph, s, c,
         fact  = sqrt(tang[0]*tang[0] + tang[1]*tang[1]),
         fact1 = sqrt((x1[0]-x2[0])*(x1[0]-x2[0]) + (x1[1]-x2[1])*(x1[1]-x2[1])),
         fact2 = sqrt((x3[0]-x2[0])*(x3[0]-x2[0]) + (x3[1]-x2[1])*(x3[1]-x2[1]));

  if (iPnt <= 0)
  {
    alph = myAcos(( (x2[0]-x1[0])*tang[0] + (x2[1]-x1[1])*tang[1]) / (fact*fact1)) * fact1/fact2;

    if (triangleArea2D(x1,x2,x3) > 0.) alph = - alph;

    c = cos(alph) / fact1 * fact * .8;
    s = sin(alph) / fact1 * fact * .8;

    tang[0] = ((x2[0]-x1[0]) * c - (x2[1]-x1[1]) * s);
    tang[1] = ((x2[0]-x1[0]) * s + (x2[1]-x1[1]) * c);

    return;
  }

  alph = myAcos(( (x3[0]-x2[0])*tang[0] + (x3[1]-x2[1])*tang[1]) / (fact*fact2)) * fact2/fact1;

  if (triangleArea2D(x1,x2,x3) < 0.) alph = - alph;

  c = cos(alph) / fact2 * fact * .8;
  s = sin(alph) / fact2 * fact * .8;

  tang[0] = ((x3[0]-x2[0]) * c - (x3[1]-x2[1]) * s);
  tang[1] = ((x3[0]-x2[0]) * s + (x3[1]-x2[1]) * c);

  return;
}





