
#include <iostream>

#include "MathBasic.h"
#include "List.h"
#include "Macro.h"
#include "Definitions.h"
#include "FunctionsProgram.h"
#include "Debug.h"
#include "MyString.h"
#include "DomainTree.h"



extern DomainTree domain;



//  Macro constructor

Macro::Macro(void) 
{
  return; 
}





Macro::Macro(int i) 
{ 
  id = i;

  if (debug) std::cout << " Macro constructor \n\n";
}





// Macro destructor

Macro::~Macro() 
{ 
  if (debug) std::cout << " Macro destructor \n\n";
}
 




//  execute macro

int Macro::exec(void)
{
  switch (id)
  {
    case   1: return macro1(*this);
    case   2: return macro2(*this);
    case   3: return macro3(*this);
    case   4: return macro4(*this);
    case   5: return macro5(*this);
    case   6: return macro6(*this);
    case   7: return macro7(*this);
    case   8: return macro8(*this);
    case   9: return macro9(*this);
    case  10: return macro10(*this);
    case  11: return macro11(*this);
    case  12: return macro12(*this);
    case  13: return macro13(*this);
    case  14: return macro14(*this);
    case  15: return macro15(*this);
    case  16: return macro16(*this);
    case  17: return macro17(*this);
    case  18: return macro18(*this);
    case  19: return macro19(*this);
    case  20: return macro20(*this);


    case 201: return macro201(*this);
    case 202: return macro202(*this);
    case 203: return macro203(*this);
    case 204: return macro204(*this);
    case 205: return macro205(*this);
    case 206: return macro206(*this);
    case 207: return macro207(*this);
  }
  return 0;
}





// overload printing to screen

std::ostream &operator<<(std::ostream &os,Macro &macro)
{
   os << "(" << macro.type << "," << macro.name << ")"; return os;
}







