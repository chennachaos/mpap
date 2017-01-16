
#include "Element2D3nodedTriangle.h"
#include "Debug.h"
#include "Plot.h"
#include "FunctionsProgram.h"
#include "FunctionsElement.h"
#include "ElementGroup.h"
#include "PropertyTypeEnum.h"
#include "MathGeom.h"
#include "Mesh.h"


extern Plot plot;


using namespace std;






Element2D3nodedTriangle::Element2D3nodedTriangle(void)
{
  ix = new int[nen()];

  if (debug) cout << " constructor Element2D3nodedTriangle\n\n";

  return;
}










Element2D3nodedTriangle::~Element2D3nodedTriangle()
{
  if (ix != NULL) delete [] ix; ix = NULL;
	
  if (debug) cout << " destructor Element2D3nodedTriangle\n\n";

  return;
}










void Element2D3nodedTriangle::plotOutline(bool defFlg)
{ 
  double &(Element::*X)(int,int);
	
  if (defFlg) X = &Element::x; else X = &Element::x0;

  for (int i=1; i<nen(); i++)
  {
    //plot.line(&((this->*X)(i,1)),&((this->*X)(i+1,1)));    
    plot.line(&(x(i,1)),&(x(i+1,1)));    
  }
  //plot.line(&((this->*X)(nen(),1)),&((this->*X)(1,1)));
  plot.line(&(x(nen(),1)),&(x(1,1)));

  return;
}










void Element2D3nodedTriangle::paint(bool defFlg)
{ 
  int n = nen();

  double &(Element::*X)(int,int);
  
  if (defFlg) X = &Element::x; else X = &Element::x0;
 
  for (int i=0; i<n; i++)
  {
    plot.polyPnt[i*3+0] = (this->*X)(i+1,1);    
    plot.polyPnt[i*3+1] = (this->*X)(i+1,2);    
  }
  plot.polyPnt[n*3+0] = (this->*X)(1,1);    
  plot.polyPnt[n*3+1] = (this->*X)(1,2);    
  
  plot.poly(n+1);
  
  return;
}










bool Element2D3nodedTriangle::forDomainType(int domType)
{
  switch (domType)
  {
    case MESH: return true;
 
    default:   return false;
  }
}










void Element2D3nodedTriangle::putLabel(char *strg, bool defFlg)
{
  int i, j;
	
  double xx[3], &(Element::*X)(int,int);

  if (defFlg) X = &Element::x; else X = &Element::x0;
  
  for (j=0; j<2; j++)  xx[j] = 0.0;
	
  for (i=0; i<3; i++)  for (j=0; j<2; j++)  xx[j] += (this->*X)(i+1,j+1); 

  for (j=0; j<2; j++) xx[j] *= 0.3333333;
  
  plot.putText(xx,strg,5);

  return;
}










void Element2D3nodedTriangle::contourPlot(int var, int indx, int nCol,
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

  plot.triangleContourPlot(&((this->*X)(1,1)),&((this->*X)(2,1)),&((this->*X)(3,1)), 
	                   (this->*U)(1,indx),(this->*U)(2,indx),(this->*U)(3,indx),
			   umn, umx, nCol); 
  return;
}










void Element2D3nodedTriangle::givePlotSequence2D(Vector<int> &plotSeq)
{
  plotSeq.free();

  for (int i=0; i<3; i++) plotSeq.append(ix[i]);

  return;
}










double Element2D3nodedTriangle::volume(bool init)
{
  if (init) 
    return  triangleArea2D(&(x0(1,1)),&(x0(2,1)),&(x0(3,1)));
	
  return  triangleArea2D(&(x(1,1)),&(x(2,1)),&(x(3,1)));
}










double Element2D3nodedTriangle::diameter(bool init)
{
  return sqrt(volume(init)) * 1.5196714;
}










void Element2D3nodedTriangle::projectVorticity(int indx)
{	
  ElementGroup *eG = (ElementGroup*) elemGrp;

  double shp[6],
	 fact = 0.5 / volume(),
	 *p   = eG->dom->p;
 
  int    i,
	 ndf_ = ndf(),
         ndm_ = ndm();
  
  if (ndf_ < ndm_) prgError(1,"Element2D3nodedTriangle::projectVorticity","ndf < ndm !");
  
  shp[0] = (x(2,2) - x(3,2)) * fact;
  shp[1] = (x(3,1) - x(2,1)) * fact;
  
  shp[2] = (x(3,2) - x(1,2)) * fact;
  shp[3] = (x(1,1) - x(3,1)) * fact;
  
  shp[4] = (x(1,2) - x(2,2)) * fact;
  shp[5] = (x(2,1) - x(1,1)) * fact;

  p[0] = 0.;
  
  for (i=0; i<3; i++)  p[0] += (u(i+1,2)*shp[i+i] - u(i+1,1)*shp[i+i+1]);

  p[1] = p[0];
  p[2] = p[0];

  return;
}








void Element2D3nodedTriangle::projectGradient(int indx, int indx2)
{
  ElementGroup *eG = (ElementGroup*) elemGrp;

  double shp[6],
	 fact = 0.5 / volume(),
	 *p   = eG->dom->p;
 
  int    i, j, k;
 
  i = min(max(1,indx),ndf());
  j = min(max(1,indx2),2); j--;
 
  shp[0] = (x(2,2) - x(3,2)) * fact;
  shp[1] = (x(3,1) - x(2,1)) * fact;
  
  shp[2] = (x(3,2) - x(1,2)) * fact;
  shp[3] = (x(1,1) - x(3,1)) * fact;
  
  shp[4] = (x(1,2) - x(2,2)) * fact;
  shp[5] = (x(2,1) - x(1,1)) * fact;

  p[0] = 0.; for (k=0; k<3; k++) p[0] += u(k+1,i) * shp[k+k+j];

  p[1] = p[0];
  p[2] = p[0];

  return;
}










void Element2D3nodedTriangle::projectNormOfGradient(int indx)
{
  ElementGroup *eG = (ElementGroup*) elemGrp;

  double shp[6],
	 fact = 0.5 / volume(),
	 *p   = eG->dom->p,
         dx = 0., dy = 0.;
 
  int    i, k;
 
  i = min(max(1,indx),ndf());
 
  shp[0] = (x(2,2) - x(3,2)) * fact;
  shp[1] = (x(3,1) - x(2,1)) * fact;
  
  shp[2] = (x(3,2) - x(1,2)) * fact;
  shp[3] = (x(1,1) - x(3,1)) * fact;
  
  shp[4] = (x(1,2) - x(2,2)) * fact;
  shp[5] = (x(2,1) - x(1,1)) * fact;

  for (k=0; k<3; k++)
  {
    dx +=  u(k+1,i) * shp[k+k];
    dy +=  u(k+1,i) * shp[k+k+1];
  }
  p[0] = sqrt(dx * dx + dy * dy);

  p[1] = p[0];
  p[2] = p[0];

  return;
}










void Element2D3nodedTriangle::getGradient(int indx, double *g)
{
  ElementGroup *eG = (ElementGroup*) elemGrp;

  double shp[6],
	 fact = 0.5 / volume();
 
  int    i, k;
 
  i = min(max(1,indx),ndf());
 
  shp[0] = (x(2,2) - x(3,2)) * fact;
  shp[1] = (x(3,1) - x(2,1)) * fact;
  
  shp[2] = (x(3,2) - x(1,2)) * fact;
  shp[3] = (x(1,1) - x(3,1)) * fact;
  
  shp[4] = (x(1,2) - x(2,2)) * fact;
  shp[5] = (x(2,1) - x(1,1)) * fact;

  g[0] = 0.;
  g[1] = 0.;

  for (k=0; k<3; k++)
  {
    g[0] +=  u(k+1,i) * shp[k+k];
    g[1] +=  u(k+1,i) * shp[k+k+1];
  }

  g[2] = sqrt(g[0] * g[0] + g[1] * g[1]);

  return;
}










void Element2D3nodedTriangle::projectNormOfGradientSquared(int indx)
{
  ElementGroup *eG = (ElementGroup*) elemGrp;

  double shp[6],
	 fact = 0.5 / volume(),
	 *p   = eG->dom->p,
         dx = 0., dy = 0.;
 
  int    i, k;
 
  i = min(max(1,indx),ndf());
 
  shp[0] = (x(2,2) - x(3,2)) * fact;
  shp[1] = (x(3,1) - x(2,1)) * fact;
  
  shp[2] = (x(3,2) - x(1,2)) * fact;
  shp[3] = (x(1,1) - x(3,1)) * fact;
  
  shp[4] = (x(1,2) - x(2,2)) * fact;
  shp[5] = (x(2,1) - x(1,1)) * fact;

  for (k=0; k<3; k++)
  {
    dx +=  u(k+1,i) * shp[k+k];
    dy +=  u(k+1,i) * shp[k+k+1];
  }
  p[0] = dx * dx + dy * dy;

  p[1] = p[0];
  p[2] = p[0];

  return;
}










int Element2D3nodedTriangle::calcStiffnessAndResidualMesh(void)
{
  ElementGroup *eG = (ElementGroup*) elemGrp;
	
  double *s      = eG->dom->s,
         *p      = eG->dom->p,
         *xl     = eG->dom->xl,
         *aleDat = &(eG->elemProp[ALETYPE]->data[0]);
	 
  int    i, j,
	 nst     = ndm()*nen(),
	 nst2    = nst  + nst,
         ndm_    = ndm(),
         nen_    = nen(),
         nale    = eG->elemProp[ALETYPE]->data.n,
	 id      = eG->elemProp[ALETYPE]->id;

  for (i=0; i<nen_; i++)
  {
    for (j=0; j<ndm_; j++) 
    {
      xl[i*ndm_+j]      = xn(i+1,j+1);
      xl[i*ndm_+j+nst]  =  d(i+1,j+1);
      xl[i*ndm_+j+nst2] = x0(i+1,j+1);
    }
  }

  // zero element residual and stiffness
  
  for (i=0; i<nst;     i++) p[i] = 0.;
  for (i=0; i<nst*nst; i++) s[i] = 0.;
 
  switch (id)
  {
    case 1: return ale2d3nodedtrianglecellcentroid_(aleDat,xl,s,p,&ndm_,&nst,&nale);
    case 2: return ale2d3nodedtriangleaspectratio_(aleDat,xl,s,p,&ndm_,&nst,&nale);
    case 3: return ale2d3nodedtrianglelinearelasticity_(aleDat,xl,s,p,&ndm_,&nst,&nale);
    case 4: return ale2d3nodedtrianglehyperelasticity_(aleDat,xl,s,p,&ndm_,&nst,&nale);
  }

  prgError(1,"Element2D3nodedTriangle::calcStiffnessAndResidualMesh","invalid ALE type!");

  return 0;
}










bool Element2D3nodedTriangle::containsPoint(double *xp, double *N)
{
  double *x0 = &(x(1,1)),
         *x1 = &(x(2,1)),
         *x2 = &(x(3,1));

  if (N == NULL) return pointInTriangle2D(x0,x1,x2,xp);

  double A0  = triangleArea2D(x0,x1,xp),
         A1  = triangleArea2D(x1,x2,xp),
         A2  = triangleArea2D(x2,x0,xp),
         dA  = 1. / (A0 + A1 + A2);
        
  N[0] = A1 * dA;
  N[1] = A2 * dA;
  N[2] = A0 * dA;

  for (int i=0; i<3; i++) if (N[i] < -0.0000001) return false;

  return true;
}



