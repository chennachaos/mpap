
#include "Element2D8nodedQuadrilateral.h"
#include "Debug.h"
#include "Plot.h"
#include "Mesh.h"
#include "MathGeom.h"


extern Plot plot;


using namespace std;








Element2D8nodedQuadrilateral::Element2D8nodedQuadrilateral(void)
{
  ix = new int[nen()];

  if (debug) cout << " constructor Element2D8nodedQuadrilateral\n\n";

  return;
}










Element2D8nodedQuadrilateral::~Element2D8nodedQuadrilateral()
{
  if (ix != NULL) delete [] ix; ix = NULL;
	
  if (debug) cout << " destructor Element2D8nodedQuadrilateral\n\n";

  return;
}










void Element2D8nodedQuadrilateral::plotOutline(bool defFlg)
{ 
  double &(Element::*X)(int,int);

  if (defFlg) X = &Element::x; else X = &Element::x0;

  for (int i=1; i<nen(); i++)
  {
    if(i<5) plot.line(&((this->*X)(i,1)),&((this->*X)(i+4,1)));  
	else plot.line(&((this->*X)(i,1)),&((this->*X)(i-3,1)));
  }
  plot.line(&((this->*X)(nen(),1)),&((this->*X)(1,1)));

  return;
}











void Element2D8nodedQuadrilateral::paint(bool defFlg)
{ 
  int n = nen();
 
  double &(Element::*X)(int,int);

  if (defFlg) X = &Element::x; else X = &Element::x0;

  for (int i=0; i<n; i++)
  {
	  if(i<4)
	  {
    plot.polyPnt[i*2*3+0] = (this->*X)(i+1,1);    
    plot.polyPnt[i*2*3+1] = (this->*X)(i+1,2);
	  }
	  else
	  {
	plot.polyPnt[(2*i-(n-1))*3+0] = (this->*X)(i+1,1);    
    plot.polyPnt[(2*i-(n-1))*3+1] = (this->*X)(i+1,2);
	  }
  }

  plot.polyPnt[n*3+0] = (this->*X)(1,1);    
  plot.polyPnt[n*3+1] = (this->*X)(1,2);    

  plot.poly(n+1);
  
  return;
}










bool Element2D8nodedQuadrilateral::forDomainType(int domType)
{
  switch (domType)
  {
    case MESH: return true;
 
    default:   return false;
  }
}










void Element2D8nodedQuadrilateral::putLabel(char *strg, bool defFlg)
{
  int i, j;
	
  double xx[3], &(Element::*X)(int,int);
  
  if (defFlg) X = &Element::x; else X = &Element::x0;

  for (j=0; j<2; j++)  xx[j] = 0.0;

  for (i=0; i<4; i++)  for (j=0; j<2; j++)  xx[j] +=  (this->*X)(i+1,j+1);

  for (j=0; j<2; j++) xx[j] *= 0.25;
  
  plot.putText(xx,strg,5);

  return;
}










void Element2D8nodedQuadrilateral::contourPlot(int var, int indx, int nCol,
		                               double umn, double umx, bool defFlg)
{
  double &(Element::*X)(int,int), &(Element::*U)(int,int), xc[3] = {0., 0., 0.}, uc = 0.;

  int i, j;
  
  if (defFlg) X = &Element::x; else X = &Element::x0;

  switch (var)
  {
    case 1: U = &Element::u; break; // degree of freedom
	    
    case 2: U = &Element::x0; break; // initial/current coordinate

    case 3: U = &Element::x; break; // current coordinate

    case 4: U = &Element::outp; break; // last projection
  }

  for (i=1; i<9; i++) 
  {
    for (j=0; j<ndm(); j++) xc[j] += (this->*X)(i,j+1);
    
    uc += (this->*U)(i,indx);
  }
  
  for (j=0; j<ndm(); j++) xc[j] *= 0.125;   uc *= 0.125;

  plot.triangleContourPlot(&((this->*X)(1,1)),&((this->*X)(5,1)),xc, 
	                   (this->*U)(1,indx),(this->*U)(5,indx),uc,
			   umn, umx, nCol); 
  
  plot.triangleContourPlot(&((this->*X)(5,1)),&((this->*X)(2,1)),xc, 
			   (this->*U)(5,indx),(this->*U)(2,indx),uc,
			   umn, umx, nCol); 
  
  plot.triangleContourPlot(&((this->*X)(2,1)),&((this->*X)(6,1)),xc, 
			   (this->*U)(2,indx),(this->*U)(6,indx),uc,
			   umn, umx, nCol); 
  
  plot.triangleContourPlot(&((this->*X)(6,1)),&((this->*X)(3,1)),xc, 
			   (this->*U)(6,indx),(this->*U)(3,indx),uc,
			   umn, umx, nCol); 

  plot.triangleContourPlot(&((this->*X)(3,1)),&((this->*X)(7,1)),xc, 
	                   (this->*U)(3,indx),(this->*U)(7,indx),uc,
			   umn, umx, nCol); 
  
  plot.triangleContourPlot(&((this->*X)(7,1)),&((this->*X)(4,1)),xc, 
			   (this->*U)(7,indx),(this->*U)(4,indx),uc,
			   umn, umx, nCol); 
  
  plot.triangleContourPlot(&((this->*X)(4,1)),&((this->*X)(8,1)),xc, 
			   (this->*U)(4,indx),(this->*U)(8,indx),uc,
			   umn, umx, nCol); 
  
  plot.triangleContourPlot(&((this->*X)(8,1)),&((this->*X)(1,1)),xc, 
			   (this->*U)(8,indx),(this->*U)(1,indx),uc,
			   umn, umx, nCol); 
  return;
}










double Element2D8nodedQuadrilateral::volume(bool init)
{
  double  &(Element::*X)(int,int), xc[2] = {0., 0.};
  int i,j;

  if (init) X = &Element::x0; else X = &Element::x;

  for (i=1; i<9; i++) 
  {
    for (j=0; j<2; j++) xc[j] += (this->*X)(i,j+1);
  }
  
  for (j=0; j<ndm(); j++) xc[j] *= 0.125;

  return  triangleArea2D(&((this->*X)(1,1)),&((this->*X)(5,1)),xc)
        + triangleArea2D(&((this->*X)(5,1)),&((this->*X)(2,1)),xc)
	+ triangleArea2D(&((this->*X)(2,1)),&((this->*X)(6,1)),xc)
	+ triangleArea2D(&((this->*X)(6,1)),&((this->*X)(3,1)),xc)
	+ triangleArea2D(&((this->*X)(3,1)),&((this->*X)(7,1)),xc)
	+ triangleArea2D(&((this->*X)(7,1)),&((this->*X)(4,1)),xc)
	+ triangleArea2D(&((this->*X)(4,1)),&((this->*X)(8,1)),xc)
	+ triangleArea2D(&((this->*X)(8,1)),&((this->*X)(1,1)),xc);
}










double Element2D8nodedQuadrilateral::diameter(bool init)
{
  return sqrt(volume(init));
}










void Element2D8nodedQuadrilateral::givePlotSequence2D(Vector<int> &plotSeq)
{
  plotSeq.free();

  for (int i=0; i<4; i++) 
  {
    plotSeq.append(ix[i]);
    plotSeq.append(ix[i+4]);
  }

  return;
}










bool Element2D8nodedQuadrilateral::containsPoint(double *xp, double *NN)
{
  double *x0   = &(x(1,1)),
         *x1   = &(x(2,1)),
         *x2   = &(x(3,1)),
         *x3   = &(x(4,1)),
         *x4   = &(x(5,1)),
         *x5   = &(x(6,1)),
         *x6   = &(x(7,1)),
         *x7   = &(x(8,1)),
         xc[2] = {.125 * (x0[0]+x1[0]+x2[0]+x3[0]+x4[0]+x5[0]+x6[0]+x7[0]),
                  .125 * (x0[1]+x1[1]+x2[1]+x3[1]+x4[1]+x5[1]+x6[1]+x7[1])},
         diag2 = ((x0[0]-x2[0])*(x0[0]-x2[0]) + (x0[1]-x2[1])*(x0[1]-x2[1])
                + (x1[0]-x3[0])*(x1[0]-x3[0]) + (x1[1]-x3[1])*(x1[1]-x3[1]));

  if (dist2(xc,xp,2) > diag2) return false;

  if (NN == NULL && (pointInTriangle2D(x4,x5,x6,xp) || pointInTriangle2D(x4,x6,x7,xp))) return true;

  double N[8], dN[8][2], xi[2] = {0.,0.}, xi2[2], R[2] = {1.,1.}, K[4], fact, tol = diag2 * 1.e-6;

  int iter = 0;

  while (iter < 10)
  {
    iter++;

    xi2[0]  = xi[0] * xi[0];
    xi2[1]  = xi[1] * xi[1];

    N[4]    = (1. - xi2[0]) * (1. -  xi[1]);
    N[5]    = (1. +  xi[0]) * (1. - xi2[1]);
    N[6]    = (1. - xi2[0]) * (1. +  xi[1]);
    N[7]    = (1. -  xi[0]) * (1. - xi2[1]);

    N[0]    = (1. - xi[0]) * (1. - xi[1]) - N[4] - N[7];
    N[1]    = (1. + xi[0]) * (1. - xi[1]) - N[5] - N[4];
    N[2]    = (1. + xi[0]) * (1. + xi[1]) - N[6] - N[5];
    N[3]    = (1. - xi[0]) * (1. + xi[1]) - N[7] - N[6];

    N[4]   += N[4];
    N[5]   += N[5];
    N[6]   += N[6];
    N[7]   += N[7];

    R[0]    =  N[0]*x0[0] + N[1]*x1[0] + N[2]*x2[0] + N[3]*x3[0]
             + N[4]*x4[0] + N[5]*x5[0] + N[6]*x6[0] + N[7]*x7[0] - 4. * xp[0];
    R[1]    =  N[0]*x0[1] + N[1]*x1[1] + N[2]*x2[1] + N[3]*x3[1]
             + N[4]*x4[1] + N[5]*x5[1] + N[6]*x6[1] + N[7]*x7[1] - 4. * xp[1];

    //cout << sqrt((R[0]*R[0]+R[1]*R[1])/diag2) << "\n";

    if (R[0]*R[0]+R[1]*R[1] < tol) 
    {
      //cout << "\n";

      if (NN != NULL) for (int i=0; i<8; i++) NN[i] = N[i] * .25;

      if (xi[0] < -1.05) return false;
      if (xi[0] > +1.05) return false;
      if (xi[1] < -1.05) return false;
      if (xi[1] > +1.05) return false;

      return true;
    }

    dN[4][0]  = (- 1. +  xi[1]) * (xi[0]+xi[0]);
    dN[4][1]  =  - 1. + xi2[0];
    dN[5][0]  =  + 1. - xi2[1];
    dN[5][1]  = (- 1. -  xi[0]) * (xi[1]+xi[1]);
    dN[6][0]  = (- 1. -  xi[1]) * (xi[0]+xi[0]);
    dN[6][1]  =  + 1. - xi2[0];
    dN[7][0]  =  - 1. + xi2[1];
    dN[7][1]  = (- 1. +  xi[0]) * (xi[1]+xi[1]);
         
    dN[0][0]  = - 1. + xi[1] - dN[4][0] - dN[7][0];
    dN[0][1]  = - 1. + xi[0] - dN[4][1] - dN[7][1];
    dN[1][0]  = + 1. - xi[1] - dN[5][0] - dN[4][0];
    dN[1][1]  = - 1. - xi[0] - dN[5][1] - dN[4][1];
    dN[2][0]  = + 1. + xi[1] - dN[6][0] - dN[5][0];
    dN[2][1]  = + 1. + xi[0] - dN[6][1] - dN[5][1];
    dN[3][0]  = - 1. - xi[1] - dN[7][0] - dN[6][0];
    dN[3][1]  = + 1. - xi[0] - dN[7][1] - dN[6][1];
         
    dN[4][0] += dN[4][0];
    dN[4][1] += dN[4][1];
    dN[5][0] += dN[5][0];
    dN[5][1] += dN[5][1];
    dN[6][0] += dN[6][0];
    dN[6][1] += dN[6][1];
    dN[7][0] += dN[7][0];
    dN[7][1] += dN[7][1];

    K[0]   =  dN[0][0]*x0[0] + dN[1][0]*x1[0] + dN[2][0]*x2[0] + dN[3][0]*x3[0]
            + dN[4][0]*x4[0] + dN[5][0]*x5[0] + dN[6][0]*x6[0] + dN[7][0]*x7[0];
    K[1]   =  dN[0][1]*x0[0] + dN[1][1]*x1[0] + dN[2][1]*x2[0] + dN[3][1]*x3[0]
            + dN[4][1]*x4[0] + dN[5][1]*x5[0] + dN[6][1]*x6[0] + dN[7][1]*x7[0];
    K[2]   =  dN[0][0]*x0[1] + dN[1][0]*x1[1] + dN[2][0]*x2[1] + dN[3][0]*x3[1]
            + dN[4][0]*x4[1] + dN[5][0]*x5[1] + dN[6][0]*x6[1] + dN[7][0]*x7[1];
    K[3]   =  dN[0][1]*x0[1] + dN[1][1]*x1[1] + dN[2][1]*x2[1] + dN[3][1]*x3[1]
            + dN[4][1]*x4[1] + dN[5][1]*x5[1] + dN[6][1]*x6[1] + dN[7][1]*x7[1];

    fact   = 1. / (K[0]*K[3]-K[1]*K[2]);

    xi[0] -= fact * (R[0] * K[3] - R[1] * K[1]);
    xi[1] -= fact * (K[0] * R[1] - K[2] * R[0]);
  }
  return false;
}










void Element2D8nodedQuadrilateral::getDistLoadFact(Vector<double> &fct, Vector<int> &nd)
{
  char fctn[] = "Element2D8nodedQuadrilateral::getDistLoadFact";

  if (nd.n != 3) prgError(1,fctn,"exactly 3 nodes required!");

  Vector<int> n, m;
  int i, l;

  for (i=0; i<3; i++)
  {
    n[i] = 0; while (n[i] < 8 && ix[n[i]] != nd[i]) n[i]++; 
    if (n[i] == 8) prgError(2,fctn,"node not found in element!");
  }

    m[0]=0; m[1]=1; m[2]=4; if (!n.containsAllOf(m))
  { m[0]=1; m[1]=2; m[2]=5; if (!n.containsAllOf(m))
  { m[0]=2; m[1]=3; m[2]=6; if (!n.containsAllOf(m))
  { m[0]=3; m[1]=0; m[2]=7; if (!n.containsAllOf(m))
    prgError(3,fctn,"wrong combination of nodes!"); }}}

  if      (n[0] > 3) { m[0] = 1; m[1] = 2; m[2] = 0; }
  else if (n[1] > 3) { m[0] = 2; m[1] = 0; m[2] = 1; }
  else               { m[0] = 0; m[1] = 1; m[2] = 2; }

  if (n[m[0]] > n[m[1]]) { i = m[0]; m[0] = m[1]; m[1] = i; }  

  //cout << nd << "\n" << n << "\n" << m << "\n\n";

  double *x1    = &(x0(n[m[0]]+1,1)),
         *x2    = &(x0(n[m[1]]+1,1)),
         *x3    = &(x0(n[m[2]]+1,1)),
         fact   = sqrt(.6),
         xi[3]  = {-fact, 0., fact},
         wgp[3] = {5./9., 8./9., 5./9.},
         N1, N2, N3, dN1, dN2, dN3, dx, dy, xi2;

  fct.free();
  for (l=0; l<3; l++)
  {
    xi2  = xi[l] * xi[l];

    N3   = 1. - xi2;
    N1   = .5 * (1. - xi[l] - N3);
    N2   = .5 * (1. + xi[l] - N3);

    dN3  = - xi[l] - xi[l];
    dN1  = .5 * ( - 1. - dN3);
    dN2  = .5 * (   1. - dN3);

    dx   = dN1*x1[0] + dN2*x2[0] + dN3*x3[0];
    dy   = dN1*x1[1] + dN2*x2[1] + dN3*x3[1];

    fact = wgp[l] * sqrt(dx*dx+dy*dy);

    fct[m[0]] += fact * N1;
    fct[m[1]] += fact * N2;
    fct[m[2]] += fact * N3;
  }

  return;
}

