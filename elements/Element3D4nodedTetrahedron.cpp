
#include "Element3D4nodedTetrahedron.h"
#include "ElementGroup.h"
#include "PropertyTypeEnum.h"
#include "FunctionsElement.h"
#include "Debug.h"
#include "Plot.h"
#include "MathBasic.h"

extern Plot plot;


using namespace std;



Element3D4nodedTetrahedron::Element3D4nodedTetrahedron(void)
{
  ix = new int[nen()];

  if (debug) cout << " constructor Element3D4nodedTetrahedron\n\n";

  return;
}




Element3D4nodedTetrahedron::~Element3D4nodedTetrahedron()
{
  if (ix != NULL) delete [] ix; ix = NULL;
	
  if (debug) cout << " destructor Element3D4nodedTetrahedron\n\n";

  return;
}





void Element3D4nodedTetrahedron::plotOutline(bool defFlg)
{ 
  MyString ch;  

  double &(Element::*X)(int,int);

  if (defFlg) X = &Element::x; else X = &Element::x0;

  plot.line(&((this->*X)(1,1)),&((this->*X)(2,1)));
  plot.line(&((this->*X)(2,1)),&((this->*X)(3,1)));
  plot.line(&((this->*X)(3,1)),&((this->*X)(1,1)));

  plot.line(&((this->*X)(1,1)),&((this->*X)(4,1)));
  plot.line(&((this->*X)(2,1)),&((this->*X)(4,1)));
  plot.line(&((this->*X)(3,1)),&((this->*X)(4,1)));


  return;
}






void Element3D4nodedTetrahedron::paint(bool defFlg)
{ 
  
  return;
}






bool Element3D4nodedTetrahedron::forDomainType(int domType)
{
  switch (domType)
  {
    case MESH: return true;
 
    default:   return false;
  }
}






void Element3D4nodedTetrahedron::putLabel(char *strg, bool defFlg)
{
/*  int i, j;
	
  double xx[3], &(Element::*X)(int,int);
  
  if (defFlg) X = &Element::x; else X = &Element::x0;

  for (j=0; j<2; j++)  xx[j] = 0.0;

  for (i=0; i<4; i++)  for (j=0; j<2; j++)  xx[j] +=  (this->*X)(i+1,j+1);

  for (j=0; j<2; j++) xx[j] *= 0.25;
  
  plot.putText(xx,strg,5);*/

  return;
}





void Element3D4nodedTetrahedron::contourPlot(int var, int indx, int nCol,
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







double Element3D4nodedTetrahedron::volume(bool init)
{
  return 1.;
}





void Element3D4nodedTetrahedron::giveFace3D(int iFace, Vector<int> &face)
{
  face.free();

  switch (iFace)
  {
    case  1: face.append(ix[0]);
             face.append(ix[2]);
             face.append(ix[1]); break;

    case  2: face.append(ix[0]);
             face.append(ix[1]);
             face.append(ix[3]); break;

    case  3: face.append(ix[3]);
             face.append(ix[2]);
             face.append(ix[0]); break;

    case  4: face.append(ix[1]);
             face.append(ix[2]);
             face.append(ix[3]); break;

    default: prgError(1,"Element3D4nodedTetrahedron::giveFace3D","invalid iFace!");
  }

  return;
}








void Element3D4nodedTetrahedron::defineBasicFace(int iFace, int i, int *fix, 
                                           unsigned int *edgeBits)
{
  fix[0] = 0;
  fix[1] = 1;
  fix[2] = 2;

  *edgeBits = 7;

  return;
}







int Element3D4nodedTetrahedron::calcStiffnessAndResidualMesh(void)
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
      xl[i*ndm_+j+nst]  = d(i+1,j+1);
      xl[i*ndm_+j+nst2] = x0(i+1,j+1);
    }
  }

  // zero element residual and stiffness
  
  for (i=0; i<nst;     i++) p[i] = 0.;
  for (i=0; i<nst*nst; i++) s[i] = 0.;




		int finite, matId; 


  switch (id)
  {
    case 2: 
							//check initial mesh 
if (aleDat[0]>1.e-14) ale3d4nodedtetrahedronaspectratiotest_(aleDat,xl,s,p,&ndm_,&nst,&nale);
  return ale3d4nodedtetrahedronaspectratio_(aleDat,xl,s,p,&ndm_,&nst,&nale);
    case 3: matId=1; finite=0;
				return ale3d4nodedtetrahedronhyperelasticity_(aleDat,xl,s,p,&matId,&finite);
    case 4: matId=2; finite=1;
				return ale3d4nodedtetrahedronhyperelasticity_(aleDat,xl,s,p,&matId,&finite);
  }

						prgError(1,"Element3D4nodedTetrahedron::calcStiffnessAndResidualMesh","invalid ALE type!");


  return 0;
}

