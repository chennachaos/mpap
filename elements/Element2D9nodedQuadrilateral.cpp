
#include "Element2D9nodedQuadrilateral.h"
#include "Debug.h"
#include "Plot.h"
#include "MathGeom.h"


extern Plot plot;


using namespace std;



Element2D9nodedQuadrilateral::Element2D9nodedQuadrilateral(void)
{
  ix = new int[nen()];

  if (debug) cout << " constructor Element2D9nodedQuadrilateral\n\n";

  return;
}




Element2D9nodedQuadrilateral::~Element2D9nodedQuadrilateral()
{
  if (ix != NULL) delete [] ix; ix = NULL;
	
  if (debug) cout << " destructor Element2D9nodedQuadrilateral\n\n";

  return;
}





void Element2D9nodedQuadrilateral::plotOutline(bool defFlg)
{ 
  double &(Element::*X)(int,int);

  if (defFlg) X = &Element::x; else X = &Element::x0;

  for (int i=1; i<nen()-1; i++)
  {
    if(i<5) plot.line(&((this->*X)(i,1)),&((this->*X)(i+4,1)));  
	else plot.line(&((this->*X)(i,1)),&((this->*X)(i-3,1)));   
  }
  plot.line(&((this->*X)(nen()-1,1)),&((this->*X)(1,1)));

  return;
}






void Element2D9nodedQuadrilateral::paint(bool defFlg)
{ 
  int n = nen()-1;
 
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






bool Element2D9nodedQuadrilateral::forDomainType(int domType)
{
  switch (domType)
  {
    case MESH: return true;
 
    default:   return false;
  }
}






void Element2D9nodedQuadrilateral::putLabel(char *strg, bool defFlg)
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

void Element2D9nodedQuadrilateral::contourPlot(int var, int indx, int nCol,
		                               double umn, double umx, bool defFlg)
{
  double &(Element::*X)(int,int), &(Element::*U)(int,int);

  int i, j;
  
  if (defFlg) X = &Element::x; else X = &Element::x0;

  switch (var)
  {
    case 1: U = &Element::u; break; // degree of freedom
	    
    case 2: U = &Element::x0; break; // initial/current coordinate

    case 3: U = &Element::x; break; // current coordinate

    case 4: U = &Element::outp; break; // last projection
  }

  plot.triangleContourPlot(&((this->*X)(1,1)),&((this->*X)(5,1)),&((this->*X)(9,1)), 
	                   (this->*U)(1,indx),(this->*U)(5,indx),(this->*U)(9,indx),
			   umn, umx, nCol); 
  
  plot.triangleContourPlot(&((this->*X)(5,1)),&((this->*X)(2,1)),&((this->*X)(9,1)), 
			   (this->*U)(5,indx),(this->*U)(2,indx),(this->*U)(9,indx),
			   umn, umx, nCol); 
  
  plot.triangleContourPlot(&((this->*X)(2,1)),&((this->*X)(6,1)),&((this->*X)(9,1)), 
			   (this->*U)(2,indx),(this->*U)(6,indx),(this->*U)(9,indx),
			   umn, umx, nCol); 
  
  plot.triangleContourPlot(&((this->*X)(6,1)),&((this->*X)(3,1)),&((this->*X)(9,1)), 
			   (this->*U)(6,indx),(this->*U)(3,indx),(this->*U)(9,indx),
			   umn, umx, nCol); 

    plot.triangleContourPlot(&((this->*X)(3,1)),&((this->*X)(7,1)),&((this->*X)(9,1)), 
	                   (this->*U)(3,indx),(this->*U)(7,indx),(this->*U)(9,indx),
			   umn, umx, nCol); 
  
  plot.triangleContourPlot(&((this->*X)(7,1)),&((this->*X)(4,1)),&((this->*X)(9,1)), 
			   (this->*U)(7,indx),(this->*U)(4,indx),(this->*U)(9,indx),
			   umn, umx, nCol); 
  
  plot.triangleContourPlot(&((this->*X)(4,1)),&((this->*X)(8,1)),&((this->*X)(9,1)), 
			   (this->*U)(4,indx),(this->*U)(8,indx),(this->*U)(9,indx),
			   umn, umx, nCol); 
  
  plot.triangleContourPlot(&((this->*X)(8,1)),&((this->*X)(1,1)),&((this->*X)(9,1)), 
			   (this->*U)(8,indx),(this->*U)(1,indx),(this->*U)(9,indx),
			   umn, umx, nCol); 

  
  return;
}







double Element2D9nodedQuadrilateral::volume(bool init)
{
  if (init)
	
    return  triangleArea2D(&(x0(1,1)),&(x0(9,1)),&(x0(5,1)))
  	  + triangleArea2D(&(x0(5,1)),&(x0(9,1)),&(x0(2,1)))
	  + triangleArea2D(&(x0(2,1)),&(x0(9,1)),&(x0(6,1)))
	  + triangleArea2D(&(x0(6,1)),&(x0(9,1)),&(x0(3,1)))
	  + triangleArea2D(&(x0(3,1)),&(x0(9,1)),&(x0(7,1)))
	  + triangleArea2D(&(x0(7,1)),&(x0(9,1)),&(x0(4,1)))
	  + triangleArea2D(&(x0(4,1)),&(x0(9,1)),&(x0(8,1)))
	  + triangleArea2D(&(x0(8,1)),&(x0(9,1)),&(x0(1,1)));
  
  return  triangleArea2D(&(x(1,1)),&(x(9,1)),&(x(5,1)))
	+ triangleArea2D(&(x(5,1)),&(x(9,1)),&(x(2,1)))
	+ triangleArea2D(&(x(2,1)),&(x(9,1)),&(x(6,1)))
	+ triangleArea2D(&(x(6,1)),&(x(9,1)),&(x(3,1)))
	+ triangleArea2D(&(x(3,1)),&(x(9,1)),&(x(7,1)))
	+ triangleArea2D(&(x(7,1)),&(x(9,1)),&(x(4,1)))
	+ triangleArea2D(&(x(4,1)),&(x(9,1)),&(x(8,1)))
	+ triangleArea2D(&(x(8,1)),&(x(9,1)),&(x(1,1)));
}








double Element2D9nodedQuadrilateral::diameter(bool init)
{
  return sqrt(volume(init));
}










void Element2D9nodedQuadrilateral::givePlotSequence2D(Vector<int> &plotSeq)
{
  plotSeq.free();

  for (int i=0; i<4; i++) 
  {
    plotSeq.append(ix[i]);
    plotSeq.append(ix[i+4]);
  }
  return;
}







bool Element2D9nodedQuadrilateral::containsPoint(double *xp, double *NN)
{
  double *x0   = &(x(1,1)),
         *x1   = &(x(2,1)),
         *x2   = &(x(3,1)),
         *x3   = &(x(4,1)),
         *x4   = &(x(5,1)),
         *x5   = &(x(6,1)),
         *x6   = &(x(7,1)),
         *x7   = &(x(8,1)),
         *x8   = &(x(9,1)),
         diag2 = ((x0[0]-x2[0])*(x0[0]-x2[0]) + (x0[1]-x2[1])*(x0[1]-x2[1])
                + (x1[0]-x3[0])*(x1[0]-x3[0]) + (x1[1]-x3[1])*(x1[1]-x3[1]));

  if (dist2(x8,xp,2) > diag2) return false;

  if (NN == NULL && (pointInTriangle2D(x4,x5,x6,xp) || pointInTriangle2D(x4,x6,x7,xp))) return true;

  double N[9], dN[9][2], xi[2] = {0.,0.}, xi2[2], R[2] = {1.,1.}, K[4], fact, tol = diag2 * 1.e-6;

  int iter = 0;

  while (iter < 10)
  {
    iter++;

    xi2[0]  = xi[0] * xi[0];
    xi2[1]  = xi[1] * xi[1];

    N[8]    = (1. - xi2[0]) * (1. - xi2[1]);

    N[4]    = (1. - xi2[0]) * (1. -  xi[1]) - N[8];
    N[5]    = (1. +  xi[0]) * (1. - xi2[1]) - N[8];
    N[6]    = (1. - xi2[0]) * (1. +  xi[1]) - N[8];
    N[7]    = (1. -  xi[0]) * (1. - xi2[1]) - N[8];

    N[0]    = (1. - xi[0]) * (1. - xi[1]) - N[4] - N[7] - N[8];
    N[1]    = (1. + xi[0]) * (1. - xi[1]) - N[5] - N[4] - N[8];
    N[2]    = (1. + xi[0]) * (1. + xi[1]) - N[6] - N[5] - N[8];
    N[3]    = (1. - xi[0]) * (1. + xi[1]) - N[7] - N[6] - N[8];

    N[4]   += N[4];
    N[5]   += N[5];
    N[6]   += N[6];
    N[7]   += N[7];

    N[8]   *= 4.;

    R[0]    =  N[0]*x0[0] + N[1]*x1[0] + N[2]*x2[0] + N[3]*x3[0] + N[4]*x4[0]
             + N[5]*x5[0] + N[6]*x6[0] + N[7]*x7[0] + N[8]*x8[0] - 4. * xp[0];
    R[1]    =  N[0]*x0[1] + N[1]*x1[1] + N[2]*x2[1] + N[3]*x3[1] + N[4]*x4[1]
             + N[5]*x5[1] + N[6]*x6[1] + N[7]*x7[1] + N[8]*x8[1] - 4. * xp[1];

    cout << sqrt((R[0]*R[0]+R[1]*R[1])/diag2) << "\n";

    if (R[0]*R[0]+R[1]*R[1] < tol) 
    {
      cout << "\n";

      if (NN != NULL) for (int i=0; i<9; i++) NN[i] = N[i] * .25;

      if (xi[0] < -1.05) return false;
      if (xi[0] > +1.05) return false;
      if (xi[1] < -1.05) return false;
      if (xi[1] > +1.05) return false;

      return true;
    }

    dN[8][0]  = (xi[0]+xi[0]) * (xi2[1] - 1.);
    dN[8][1]  = (xi[1]+xi[1]) * (xi2[0] - 1.);

    dN[4][0]  = (- 1. +  xi[1]) * (xi[0]+xi[0]) - dN[8][0];
    dN[4][1]  =  - 1. + xi2[0]                  - dN[8][1];
    dN[5][0]  =  + 1. - xi2[1]                  - dN[8][0];
    dN[5][1]  = (- 1. -  xi[0]) * (xi[1]+xi[1]) - dN[8][1];
    dN[6][0]  = (- 1. -  xi[1]) * (xi[0]+xi[0]) - dN[8][0];
    dN[6][1]  =  + 1. - xi2[0]                  - dN[8][1];
    dN[7][0]  =  - 1. + xi2[1]                  - dN[8][0];
    dN[7][1]  = (- 1. +  xi[0]) * (xi[1]+xi[1]) - dN[8][1];
         
    dN[0][0]  = - 1. + xi[1] - dN[4][0] - dN[7][0] - dN[8][0];
    dN[0][1]  = - 1. + xi[0] - dN[4][1] - dN[7][1] - dN[8][1];
    dN[1][0]  = + 1. - xi[1] - dN[5][0] - dN[4][0] - dN[8][0];
    dN[1][1]  = - 1. - xi[0] - dN[5][1] - dN[4][1] - dN[8][1];
    dN[2][0]  = + 1. + xi[1] - dN[6][0] - dN[5][0] - dN[8][0];
    dN[2][1]  = + 1. + xi[0] - dN[6][1] - dN[5][1] - dN[8][1];
    dN[3][0]  = - 1. - xi[1] - dN[7][0] - dN[6][0] - dN[8][0];
    dN[3][1]  = + 1. - xi[0] - dN[7][1] - dN[6][1] - dN[8][1];
         
    dN[4][0] += dN[4][0];
    dN[4][1] += dN[4][1];
    dN[5][0] += dN[5][0];
    dN[5][1] += dN[5][1];
    dN[6][0] += dN[6][0];
    dN[6][1] += dN[6][1];
    dN[7][0] += dN[7][0];
    dN[7][1] += dN[7][1];

    dN[8][0] *= 4.;
    dN[8][1] *= 4.;

    K[0]   =  dN[0][0]*x0[0] + dN[1][0]*x1[0] + dN[2][0]*x2[0] + dN[3][0]*x3[0]
            + dN[4][0]*x4[0] + dN[5][0]*x5[0] + dN[6][0]*x6[0] + dN[7][0]*x7[0] + dN[8][0]*x8[0];
    K[1]   =  dN[0][1]*x0[0] + dN[1][1]*x1[0] + dN[2][1]*x2[0] + dN[3][1]*x3[0]
            + dN[4][1]*x4[0] + dN[5][1]*x5[0] + dN[6][1]*x6[0] + dN[7][1]*x7[0] + dN[8][1]*x8[0];
    K[2]   =  dN[0][0]*x0[1] + dN[1][0]*x1[1] + dN[2][0]*x2[1] + dN[3][0]*x3[1]
            + dN[4][0]*x4[1] + dN[5][0]*x5[1] + dN[6][0]*x6[1] + dN[7][0]*x7[1] + dN[8][0]*x8[1];
    K[3]   =  dN[0][1]*x0[1] + dN[1][1]*x1[1] + dN[2][1]*x2[1] + dN[3][1]*x3[1]
            + dN[4][1]*x4[1] + dN[5][1]*x5[1] + dN[6][1]*x6[1] + dN[7][1]*x7[1] + dN[8][1]*x8[1];

    fact   = 1. / (K[0]*K[3]-K[1]*K[2]);

    xi[0] -= fact * (R[0] * K[3] - R[1] * K[1]);
    xi[1] -= fact * (K[0] * R[1] - K[2] * R[0]);
  }
  return false;
}

