
#include "ElementGroup.h"
#include "Debug.h"





ElementGroup::ElementGroup(void)
{
  return;
}






ElementGroup::ElementGroup(Domain* domPtr)
{
  dom = domPtr;

  ndf = 0;
  ndm = 0;
  nen = 0;
  
  nivGP = 0;
  nGP   = 0;
  
  elemProp  = NULL;
  nElemProp = 0;

  closed = true;
  
  if (debug) std::cout << " constructor ElementGroup\n\n";
  
  return;
}




ElementGroup::~ElementGroup()
{
  if (elemProp != NULL) delete [] elemProp;  elemProp = NULL;
  
  if (debug) std::cout << " destructor ElementGroup\n\n";
  
  return;
}




void ElementGroup::addElemProp(int type, PropertyItem *eP)
{
  if (type < nElemProp)  {  elemProp[type] = eP;   return;  }

  PropertyItem **tmp = new PropertyItem* [type+1];
  
  int i;
  
  for (i=0; i<nElemProp; i++)    tmp[i] = elemProp[i];

  for (i=nElemProp; i<type; i++) tmp[i] = NULL;

  tmp[i] = eP;

  delete [] elemProp; elemProp = tmp;  

  nElemProp = type + 1;
	
  return;
}


