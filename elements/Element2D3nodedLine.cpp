
#include "Element2D3nodedLine.h"
#include "Debug.h"
#include "Plot.h"


extern Plot plot;


using namespace std;



Element2D3nodedLine::Element2D3nodedLine(void)
{
  ix = new int[nen()];
  
  if (debug) cout << " constructor Element2D3nodedLine\n\n";

  return;
}




Element2D3nodedLine::~Element2D3nodedLine()
{
  if (ix != NULL) delete [] ix; ix = NULL;
	
  if (debug) cout << " destructor Element2D3nodedLine\n\n";

  return;
}





void Element2D3nodedLine::plotOutline(bool defFlg)
{ 
  double &(Element::*X)(int,int);
	
  if (defFlg) X = &Element::x; else X = &Element::x0;

  plot.line(&((this->*X)(1,1)),&((this->*X)(2,1)));
  plot.line(&((this->*X)(2,1)),&((this->*X)(3,1)));

  return;
}






void Element2D3nodedLine::paint(bool defFlg)
{ 
  plotOutline(defFlg);
  
  return;
}






bool Element2D3nodedLine::forDomainType(int domType)
{
  switch (domType)
  {
    case MESH: return true;
 
    default:   return false;
  }
}





void Element2D3nodedLine::putLabel(char *strg, bool defFlg)
{
  int i, j;
	
  double &(Element::*X)(int,int);

  if (defFlg) X = &Element::x; else X = &Element::x0;
  
  plot.putText(&((this->*X)(2,1)),strg,5);

  return;
}






void Element2D3nodedLine::givePlotSequence2D(Vector<int> &plotSeq)
{
  plotSeq.free();
      
  plotSeq[0] = ix[0];
  plotSeq[1] = ix[1];
  plotSeq[2] = ix[2];

  return;
}

