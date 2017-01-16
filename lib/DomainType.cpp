
#include <iostream>

#include "Definitions.h"
#include "FunctionsProgram.h"
#include "DomainType.h"



DomainType::DomainType(int t)
{
  char *domType [] = DOMAIN_TYPE_NAMES;
      
  tid  = t;

  type = domType[t];

  return;
}





DomainType::~DomainType(void)
{
  type.free();
  
  dom.free();
  
  key.free();
  var.free();
  fct.free();
  
  return;
}


 


char *DomainType::printName(void)
{  
  return type.asCharArray();  
}





bool DomainType::operator==(const DomainType &item)
{
  if (tid == item.tid) return true;

  return false;
}




