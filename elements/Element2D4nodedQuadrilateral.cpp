
#include "Element2D4nodedQuadrilateral.h"
#include "Debug.h"
#include "Plot.h"
#include "MathGeom.h"
#include "FunctionsElement.h"
#include "ElementGroup.h"
#include "PropertyTypeEnum.h"
#include "Mesh.h"



extern Plot plot;


using namespace std;






Element2D4nodedQuadrilateral::Element2D4nodedQuadrilateral(void)
{
  ix = new int[nen()];

  if (debug) cout << " constructor Element2D4nodedQuadrilateral\n\n";

  return;
}










Element2D4nodedQuadrilateral::~Element2D4nodedQuadrilateral()
{
  if (ix != NULL) delete [] ix; ix = NULL;
	
  if (debug) cout << " destructor Element2D4nodedQuadrilateral\n\n";

  return;
}










void Element2D4nodedQuadrilateral::plotOutline(bool defFlg)
{ 
  double &(Element::*X)(int,int);

  if (defFlg) X = &Element::x; else X = &Element::x0;

  for (int i=1; i<nen(); i++)
  {
    plot.line(&((this->*X)(i,1)),&((this->*X)(i+1,1)));    
  }
  plot.line(&((this->*X)(nen(),1)),&((this->*X)(1,1)));

  return;
}










void Element2D4nodedQuadrilateral::paint(bool defFlg)
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










bool Element2D4nodedQuadrilateral::forDomainType(int domType)
{
  switch (domType)
  {
    case MESH: return true;
 
    default:   return false;
  }
}










void Element2D4nodedQuadrilateral::putLabel(char *strg, bool defFlg)
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










void Element2D4nodedQuadrilateral::contourPlot(int var, int indx, int nCol,
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

  for (i=1; i<5; i++) 
  {
    for (j=0; j<ndm(); j++) xc[j] += (this->*X)(i,j+1);
    
    uc += (this->*U)(i,indx);
  }
  
  for (j=0; j<ndm(); j++) xc[j] *= 0.25;   uc *= 0.25;

  plot.triangleContourPlot(&((this->*X)(1,1)),&((this->*X)(2,1)),xc, 
	                   (this->*U)(1,indx),(this->*U)(2,indx),uc,
			   umn, umx, nCol); 
  
  plot.triangleContourPlot(&((this->*X)(2,1)),&((this->*X)(3,1)),xc, 
			   (this->*U)(2,indx),(this->*U)(3,indx),uc,
			   umn, umx, nCol); 
  
  plot.triangleContourPlot(&((this->*X)(3,1)),&((this->*X)(4,1)),xc, 
			   (this->*U)(3,indx),(this->*U)(4,indx),uc,
			   umn, umx, nCol); 
  
  plot.triangleContourPlot(&((this->*X)(4,1)),&((this->*X)(1,1)),xc, 
			   (this->*U)(4,indx),(this->*U)(1,indx),uc,
			   umn, umx, nCol); 

  //MyString ch; ch.inputKeepIfReturn();
  
  return;
}










double Element2D4nodedQuadrilateral::volume(bool init)
{
  if (init) 
    return  triangleArea2D(&(x0(1,1)),&(x0(2,1)),&(x0(3,1)))
          + triangleArea2D(&(x0(1,1)),&(x0(3,1)),&(x0(4,1)));
	
  return  triangleArea2D(&(x(1,1)),&(x(2,1)),&(x(3,1)))
        + triangleArea2D(&(x(1,1)),&(x(3,1)),&(x(4,1)));
}










double Element2D4nodedQuadrilateral::diameter(bool init)
{
  return sqrt(volume(init));
}










void Element2D4nodedQuadrilateral::givePlotSequence2D(Vector<int> &plotSeq)
{
  plotSeq.free();

  for (int i=0; i<4; i++) plotSeq.append(ix[i]);

  return;
}










int Element2D4nodedQuadrilateral::calcStiffnessAndResidualMesh(void)
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
    case 2: return ale2d4nodedquadrilateralaspectratio_(aleDat,xl,s,p,&ndm_,&nst,&nale);
  }

  prgError(1,"Element2D4nodedQuadrilateral::calcStiffnessAndResidualMesh","invalid ALE type!");

  return 0;
}










bool Element2D4nodedQuadrilateral::containsPoint(double *xp, double *NN)
{
  double *x0   = &(x(1,1)),
         *x1   = &(x(2,1)),
         *x2   = &(x(3,1)),
         *x3   = &(x(4,1));

  if (NN == NULL) return (pointInTriangle2D(x0,x1,x2,xp) || pointInTriangle2D(x0,x2,x3,xp));

  double N[4], dN[4][2], xi[2] = {0.,0.}, R[2] = {1.,1.}, K[4], fact, tol = diameter();

  tol *= (tol * 1.e-6);

  int iter = 0;

  while (iter < 10)
  {
    iter++;

    N[0]   = (1. - xi[0]) * (1. - xi[1]);
    N[1]   = (1. + xi[0]) * (1. - xi[1]);
    N[2]   = (1. + xi[0]) * (1. + xi[1]);
    N[3]   = (1. - xi[0]) * (1. + xi[1]);

    R[0]   = N[0]*x0[0] + N[1]*x1[0] + N[2]*x2[0] + N[3]*x3[0] - 4.*xp[0];
    R[1]   = N[0]*x0[1] + N[1]*x1[1] + N[2]*x2[1] + N[3]*x3[1] - 4.*xp[1];

    //cout << sqrt(R[0]*R[0]+R[1]*R[1])/diameter() << "\n";

    if (R[0]*R[0]+R[1]*R[1] < tol) 
    {
      //cout << "\n";

      for (int i=0; i<4; i++) NN[i] = N[i] * .25;

      if (xi[0] < -1.05) return false;
      if (xi[0] > +1.05) return false;
      if (xi[1] < -1.05) return false;
      if (xi[1] > +1.05) return false;

      return true;
    }

    dN[0][0] = - 1. + xi[1];
    dN[0][1] = - 1. + xi[0];
    dN[1][0] = + 1. - xi[1];
    dN[1][1] = - 1. - xi[0];
    dN[2][0] = + 1. + xi[1];
    dN[2][1] = + 1. + xi[0];
    dN[3][0] = - 1. - xi[1];
    dN[3][1] = + 1. - xi[0];

    K[0]   = dN[0][0]*x0[0] + dN[1][0]*x1[0] + dN[2][0]*x2[0] + dN[3][0]*x3[0];
    K[1]   = dN[0][1]*x0[0] + dN[1][1]*x1[0] + dN[2][1]*x2[0] + dN[3][1]*x3[0];
    K[2]   = dN[0][0]*x0[1] + dN[1][0]*x1[1] + dN[2][0]*x2[1] + dN[3][0]*x3[1];
    K[3]   = dN[0][1]*x0[1] + dN[1][1]*x1[1] + dN[2][1]*x2[1] + dN[3][1]*x3[1];

    fact   = 1. / (K[0]*K[3]-K[1]*K[2]);

    xi[0] -= fact * (R[0] * K[3] - R[1] * K[1]);
    xi[1] -= fact * (K[0] * R[1] - K[2] * R[0]);
  }

  return false;
}










void Element2D4nodedQuadrilateral::getDistLoadFact(Vector<double> &fct, Vector<int> &nd)
{


  return;
}





