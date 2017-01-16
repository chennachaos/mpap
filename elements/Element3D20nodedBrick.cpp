
#include "Element3D20nodedBrick.h"
#include "Debug.h"
#include "Plot.h"
#include "MathBasic.h"

extern Plot plot;


using namespace std;



Element3D20nodedBrick::Element3D20nodedBrick(void)
{
  ix = new int[nen()];

  if (debug) cout << " constructor Element3D20nodedBrick\n\n";

  return;
}




Element3D20nodedBrick::~Element3D20nodedBrick()
{
  if (ix != NULL) delete [] ix; ix = NULL;
	
  if (debug) cout << " destructor Element3D20nodedBrick\n\n";

  return;
}





void Element3D20nodedBrick::plotOutline(bool defFlg)
{ 

  MyString ch; 

  int i;

  double &(Element::*X)(int,int);

  if (defFlg) X = &Element::x; else X = &Element::x0;

 	  for (i=0; i<10; i++)  
	  {
	  plot.line(&((this->*X)(i+1,1)),&((this->*X)(i+9,1)));

	  if(i<3) 
	  {
	  plot.line(&((this->*X)(i+9,1)),&((this->*X)(i+2,1)));
	  plot.line(&((this->*X)(i+13,1)),&((this->*X)(i+6,1)));
	  }
	  if(i<4) 
	  {
	  plot.line(&((this->*X)(i+1,1)),&((this->*X)(i+17,1)));
	  plot.line(&((this->*X)(i+17,1)),&((this->*X)(i+5,1)));
	  }
	  }
   plot.line(&((this->*X)(16,1)),&((this->*X)(5,1)));
   plot.line(&((this->*X)(12,1)),&((this->*X)(1,1)));
  

  return;
}






void Element3D20nodedBrick::paint(bool defFlg)
{ 
  
  return;
}






bool Element3D20nodedBrick::forDomainType(int domType)
{
  switch (domType)
  {
    case MESH: return true;
 
    default:   return false;
  }
}






void Element3D20nodedBrick::putLabel(char *strg, bool defFlg)
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





void Element3D20nodedBrick::contourPlot(int var, int indx, int nCol,
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







double Element3D20nodedBrick::volume(bool init)
{
  return 1.;

}





void Element3D20nodedBrick::giveFace3D(int iFace, Vector<int> &face)
{
  face.free();

  switch (iFace)
  {
    case  1: face.append(ix[0]);
             face.append(ix[8]);
             face.append(ix[1]);
             face.append(ix[17]);
             face.append(ix[5]);
             face.append(ix[12]);
             face.append(ix[4]); 
             face.append(ix[16]); break;

    case  2: face.append(ix[1]);
             face.append(ix[9]);
             face.append(ix[2]);
             face.append(ix[18]);
             face.append(ix[6]);
             face.append(ix[13]);
             face.append(ix[5]); 			 
             face.append(ix[17]); break;

    case  3: face.append(ix[2]);
             face.append(ix[10]);
             face.append(ix[3]);
             face.append(ix[19]);
             face.append(ix[7]);
             face.append(ix[14]);
             face.append(ix[6]); 			 
             face.append(ix[18]); break;

    case  4: face.append(ix[3]);
             face.append(ix[11]);
             face.append(ix[0]);
             face.append(ix[16]);
             face.append(ix[4]);
             face.append(ix[15]);
             face.append(ix[7]); 			
             face.append(ix[19]); break;

    case  5: face.append(ix[3]);
             face.append(ix[10]);
             face.append(ix[2]);
             face.append(ix[9]);
             face.append(ix[1]);
             face.append(ix[8]);
             face.append(ix[0]); 			 
             face.append(ix[11]); break;

    case  6: face.append(ix[4]);
             face.append(ix[12]);
             face.append(ix[5]);
             face.append(ix[13]);
             face.append(ix[6]);
             face.append(ix[14]);
             face.append(ix[7]); 			 
             face.append(ix[15]); break;

    default: prgError(1,"Element3D20nodedBrick::giveFace3D","invalid iFace!");
  }

  return;
}








void Element3D20nodedBrick::defineBasicFace(int iFace, int i, int *fix, 
                                           unsigned int *edgeBits)
{
  switch (i)
  {
   case 5:  fix[0] = 7;
            fix[1] = 0;
            fix[2] = 1;
           
            *edgeBits = 3;
           
            break;

   case 4:  fix[0] = 1;
            fix[1] = 2;
            fix[2] = 3;
           
            *edgeBits = 3;
           
            break;

   case 3:  fix[0] = 3;
            fix[1] = 4;
            fix[2] = 5;
           
            *edgeBits = 3;
           
            break;

   case 2:  fix[0] = 5;
            fix[1] = 6;
            fix[2] = 7;
           
            *edgeBits = 3;
           
            break;

   case 1:  fix[0] = 1;
            fix[1] = 3;
            fix[2] = 5;
           
            *edgeBits = 0;
           
            break;

   case 0:  fix[0] = 1;
            fix[1] = 5;
            fix[2] = 7;

            *edgeBits = 0;

            break;
  }

  return;
}




