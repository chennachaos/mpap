
#include "Element2D6nodedTriangle.h"
#include "Debug.h"
#include "Plot.h"
#include "MathGeom.h"


extern Plot plot;


using namespace std;



Element2D6nodedTriangle::Element2D6nodedTriangle(void)
{
  ix = new int[nen()];

  if (debug) cout << " constructor Element2D6nodedTriangle\n\n";

  return;
}




Element2D6nodedTriangle::~Element2D6nodedTriangle()
{
  if (ix != NULL) delete [] ix; ix = NULL;
	
  if (debug) cout << " destructor Element2D6nodedTriangle\n\n";

  return;
}





void Element2D6nodedTriangle::plotOutline(bool defFlg)
{ 
  double &(Element::*X)(int,int);

  if (defFlg) X = &Element::x; else X = &Element::x0;

  for (int i=1; i<4; i++)
  {
     plot.line(&((this->*X)(i,1)),&((this->*X)(i+3,1)));  
	if(i<3) plot.line(&((this->*X)(i+3,1)),&((this->*X)(i+1,1)));  

  }
  plot.line(&((this->*X)(nen(),1)),&((this->*X)(1,1)));

  return;
}






void Element2D6nodedTriangle::paint(bool defFlg)
{ 
  int n = nen();
 
  double &(Element::*X)(int,int);

  if (defFlg) X = &Element::x; else X = &Element::x0;

  for (int i=0; i<n; i++)
  {
	  if(i<3)
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






bool Element2D6nodedTriangle::forDomainType(int domType)
{
  switch (domType)
  {
    case MESH: return true;
 
    default:   return false;
  }
}






void Element2D6nodedTriangle::putLabel(char *strg, bool defFlg)
{
  int i, j;
	
  double xx[3], &(Element::*X)(int,int);
  
  if (defFlg) X = &Element::x; else X = &Element::x0;

  for (j=0; j<2; j++)  xx[j] = 0.0;

  for (i=0; i<3; i++)  for (j=0; j<2; j++)  xx[j] +=  (this->*X)(i+1,j+1);

  for (j=0; j<2; j++) xx[j] /= 3.0;
  
  plot.putText(xx,strg,5);

  return;
}


void Element2D6nodedTriangle::contourPlot(int var, int indx, int nCol,
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

  for (i=1; i<7; i++) 
  {
    for (j=0; j<ndm(); j++) xc[j] += (this->*X)(i,j+1);
    
    uc += (this->*U)(i,indx);
  }
  
  for (j=0; j<ndm(); j++) xc[j] /= 6.0;   uc /= 6.0;

  plot.triangleContourPlot(&((this->*X)(1,1)),&((this->*X)(4,1)),xc, 
	                   (this->*U)(1,indx),(this->*U)(4,indx),uc,
			   umn, umx, nCol); 
  
  plot.triangleContourPlot(&((this->*X)(4,1)),&((this->*X)(2,1)),xc, 
			   (this->*U)(4,indx),(this->*U)(2,indx),uc,
			   umn, umx, nCol); 
  
  plot.triangleContourPlot(&((this->*X)(2,1)),&((this->*X)(5,1)),xc, 
			   (this->*U)(2,indx),(this->*U)(5,indx),uc,
			   umn, umx, nCol); 
  
  plot.triangleContourPlot(&((this->*X)(5,1)),&((this->*X)(3,1)),xc, 
			   (this->*U)(5,indx),(this->*U)(3,indx),uc,
			   umn, umx, nCol); 

    plot.triangleContourPlot(&((this->*X)(3,1)),&((this->*X)(6,1)),xc, 
	                   (this->*U)(3,indx),(this->*U)(6,indx),uc,
			   umn, umx, nCol); 
  
   plot.triangleContourPlot(&((this->*X)(6,1)),&((this->*X)(1,1)),xc, 
			   (this->*U)(6,indx),(this->*U)(1,indx),uc,
			   umn, umx, nCol); 
  
  return;
}







double Element2D6nodedTriangle::volume(bool init)
{
  double  &(Element::*X)(int,int), xc[2] = {0., 0.};
  int i,j;

  if (init) X = &Element::x0; else X = &Element::x;

  for (i=1; i<7; i++) 
  {
    for (j=0; j<2; j++) xc[j] += (this->*X)(i,j+1);
  }
  
  for (j=0; j<ndm(); j++) xc[j] /= 6.0;

  return  triangleArea2D(&((this->*X)(1,1)),&((this->*X)(4,1)),xc)
        + triangleArea2D(&((this->*X)(4,1)),&((this->*X)(2,1)),xc)
	+ triangleArea2D(&((this->*X)(2,1)),&((this->*X)(5,1)),xc)
	+ triangleArea2D(&((this->*X)(5,1)),&((this->*X)(3,1)),xc)
	+ triangleArea2D(&((this->*X)(3,1)),&((this->*X)(6,1)),xc)
	+ triangleArea2D(&((this->*X)(6,1)),&((this->*X)(1,1)),xc);
}







double Element2D6nodedTriangle::diameter(bool init)
{
  return sqrt(volume(init)) * 1.5196714;
}










void Element2D6nodedTriangle::givePlotSequence2D(Vector<int> &plotSeq)
{
  plotSeq.free();

  for (int i=0; i<3; i++) 
  {
    plotSeq.append(ix[i]);
    plotSeq.append(ix[i+3]);
  }

  return;
}












bool Element2D6nodedTriangle::containsPoint(double *xp, double *NN)
{
  double *x0   = &(x(1,1)),
         *x1   = &(x(2,1)),
         *x2   = &(x(3,1)),
         *x3   = &(x(4,1)),
         *x4   = &(x(5,1)),
         *x5   = &(x(6,1)),
         xc[2] = {.1666667 * (x0[0]+x1[0]+x2[0]+x3[0]+x4[0]+x5[0]),
                  .1666667 * (x0[1]+x1[1]+x2[1]+x3[1]+x4[1]+x5[1])},
         diag2 = 0.6 * (dist2(xc,x0,2) + dist2(xc,x1,2) + dist2(xc,x2,2));

  if (dist2(xc,xp,2) > diag2) return false;

  if (NN == NULL && pointInTriangle2D(x3,x4,x5,xp)) return true;

  double N[6], dN[6][2], xi[2] = {0.,0.}, xi2[2], R[2] = {1.,1.}, K[4], fact, fact2, 
         tol = diag2 * 1.e-8;

  int iter = 0;

  while (iter < 10)
  {
    iter++;

    xi2[0]  = xi[0] * xi[0];
    xi2[1]  = xi[1] * xi[1];

    fact    = 1. - xi[0] - xi[1];
    fact2   = fact * fact; 

    N[0]    = xi2[0] + xi2[0] - xi[0];
    N[1]    = xi2[1] + xi2[1] - xi[1];
    N[2]    = fact2  + fact2  - fact;
    N[3]    = 4. * xi[0] * xi[1];
    N[4]    = 4. * xi[1] * fact;
    N[5]    = 4. * xi[0] * fact;

    R[0]    = N[0]*x0[0] + N[1]*x1[0] + N[2]*x2[0] + N[3]*x3[0] + N[4]*x4[0] + N[5]*x5[0] - xp[0];
    R[1]    = N[0]*x0[1] + N[1]*x1[1] + N[2]*x2[1] + N[3]*x3[1] + N[4]*x4[1] + N[5]*x5[1] - xp[1];

    //cout << sqrt((R[0]*R[0]+R[1]*R[1])/diag2) << "\n";

    if (R[0]*R[0]+R[1]*R[1] < tol) 
    {
      //cout << "\n";

      if (NN != NULL) for (int i=0; i<6; i++) NN[i] = N[i];

      if (xi[0]       < -0.02) return false;
      if (xi[0]+xi[1] > +1.02) return false;
      if (xi[1]       < -0.02) return false;

      return true;
    }

    fact2     = - fact - fact;

    dN[0][0]  = 4. * xi[0] - 1.;
    dN[0][1]  = 0.;
    dN[1][0]  = 0.;
    dN[1][1]  = 4. * xi[1] - 1.;
    dN[2][0]  = fact2 + fact2 + 1.;
    dN[2][1]  = fact2 + fact2 + 1.;
    dN[3][0]  = 4. * xi[1]; 
    dN[3][1]  = 4. * xi[0];
    dN[4][0]  = - 4. * xi[1];
    dN[4][1]  = 4. * (fact - xi[1]);
    dN[5][0]  = 4. * (fact - xi[0]);
    dN[5][1]  = - 4. * xi[0];
         
    K[0]   =  dN[0][0]*x0[0] + dN[1][0]*x1[0] + dN[2][0]*x2[0] + dN[3][0]*x3[0]
            + dN[4][0]*x4[0] + dN[5][0]*x5[0];
    K[1]   =  dN[0][1]*x0[0] + dN[1][1]*x1[0] + dN[2][1]*x2[0] + dN[3][1]*x3[0]
            + dN[4][1]*x4[0] + dN[5][1]*x5[0];
    K[2]   =  dN[0][0]*x0[1] + dN[1][0]*x1[1] + dN[2][0]*x2[1] + dN[3][0]*x3[1]
            + dN[4][0]*x4[1] + dN[5][0]*x5[1];
    K[3]   =  dN[0][1]*x0[1] + dN[1][1]*x1[1] + dN[2][1]*x2[1] + dN[3][1]*x3[1]
            + dN[4][1]*x4[1] + dN[5][1]*x5[1];

    fact   = 1. / (K[0]*K[3]-K[1]*K[2]);

    xi[0] -= fact * (R[0] * K[3] - R[1] * K[1]);
    xi[1] -= fact * (K[0] * R[1] - K[2] * R[0]);
  }
  return false;
}

