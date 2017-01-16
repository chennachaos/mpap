
#include "Element1D2nodedLine.h"
#include "Debug.h"
#include "Plot.h"


extern Plot plot;


using namespace std;



Element1D2nodedLine::Element1D2nodedLine(void)
{
  ix = new int[nen()];

  if (debug) cout << " constructor Element1D2nodedLine\n\n";

  return;
}




Element1D2nodedLine::~Element1D2nodedLine()
{
  if (ix != NULL) delete [] ix; ix = NULL;
	
  if (debug) cout << " destructor Element1D2nodedLine\n\n";

  return;
}





void Element1D2nodedLine::plotOutline(bool defFlg)
{ 
//  double &(Element::*X)(int,int);
	
//  X = &Element::x0;

//  plot.line(&((this->*X)(nen(),1)),&((this->*X)(1,1)));

  return;
}






void Element1D2nodedLine::paint(bool defFlg)
{ 
  plotOutline(defFlg);
  
  return;
}






bool Element1D2nodedLine::forDomainType(int domType)
{
  switch (domType)
  {
    case MESH: return true;
 
    default:   return false;
  }
}





void Element1D2nodedLine::putLabel(char *strg, bool defFlg)
{
//  int i, j;
	
//  double xx[3], &(Element::*X)(int,int);

//  if (defFlg) X = &Element::x; else X = &Element::x0;
  
//  for (j=0; j<2; j++)  xx[j] = 0.0;
	
//  for (i=0; i<2; i++)  for (j=0; j<2; j++)  xx[j] += (this->*X)(i+1,j+1); 

//  for (j=0; j<2; j++) xx[j] *= 0.5;
  
//  plot.putText(xx,strg,5);

  return;
}







