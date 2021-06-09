
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
    case 208: return macro208(*this); 
    case 209: return macro209(*this); 
    case 210: return macro210(*this); 
    case 211: return macro211(*this); 
    case 212: return macro212(*this); 
    case 213: return macro213(*this); 
    case 214: return macro214(*this); 
    case 215: return macro215(*this); 
    case 216: return macro216(*this); 
    case 217: return macro217(*this); 
    case 218: return macro218(*this); 
    case 219: return macro219(*this); 
    case 220: return macro220(*this); 
    case 221: return macro221(*this); 
    case 222: return macro222(*this); 
    case 223: return macro223(*this); 
    case 224: return macro224(*this); 
    case 225: return macro225(*this); 
    case 226: return macro226(*this); 
    case 227: return macro227(*this); 
    case 228: return macro228(*this); 
    case 229: return macro229(*this); 
    case 230: return macro230(*this); 
  }
  return 0;
}





// overload printing to screen

std::ostream &operator<<(std::ostream &os,Macro &macro)
{
   os << "(" << macro.type << "," << macro.name << ")"; return os;
}







