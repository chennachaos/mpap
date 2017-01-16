
#include "Element2D2nodedLine.h"
#include "Debug.h"
#include "Plot.h"


extern Plot plot;


using namespace std;



Element2D2nodedLine::Element2D2nodedLine(void)
{
  ix = new int[nen()];
  
  if (debug) cout << " constructor Element2D2nodedLine\n\n";

  return;
}




Element2D2nodedLine::~Element2D2nodedLine()
{
  if (ix != NULL) delete [] ix; ix = NULL;
	
  if (debug) cout << " destructor Element2D2nodedLine\n\n";

  return;
}





void Element2D2nodedLine::plotOutline(bool defFlg)
{ 
  double &(Element::*X)(int,int);
	
  if (defFlg) X = &Element::x; else X = &Element::x0;

  plot.line(&((this->*X)(nen(),1)),&((this->*X)(1,1)));

  return;
}






void Element2D2nodedLine::paint(bool defFlg)
{ 
  plotOutline(defFlg);
  
  return;
}






bool Element2D2nodedLine::forDomainType(int domType)
{
  switch (domType)
  {
    case MESH: return true;
 
    default:   return false;
  }
}





void Element2D2nodedLine::putLabel(char *strg, bool defFlg)
{
  int i, j;
	
  double xx[3], &(Element::*X)(int,int);

  if (defFlg) X = &Element::x; else X = &Element::x0;
  
  for (j=0; j<2; j++)  xx[j] = 0.0;
	
  for (i=0; i<2; i++)  for (j=0; j<2; j++)  xx[j] += (this->*X)(i+1,j+1); 

  for (j=0; j<2; j++) xx[j] *= 0.5;
  
  plot.putText(xx,strg,5);

  return;
}






void Element2D2nodedLine::givePlotSequence2D(Vector<int> &plotSeq)
{
  plotSeq.free();
      
  plotSeq.append(ix[0]);
  plotSeq.append(ix[1]);

  return;
}

