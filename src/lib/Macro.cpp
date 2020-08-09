
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
    case  21: return macro21(*this); 
    case  22: return macro22(*this); 
    case  23: return macro23(*this); 
    case  24: return macro24(*this); 
    case  25: return macro25(*this); 
    case  26: return macro26(*this); 
    case  27: return macro27(*this); 
    case  28: return macro28(*this); 
    case  29: return macro29(*this); 
    case  30: return macro30(*this); 
    case  31: return macro31(*this); 
    case  32: return macro32(*this); 
    case  33: return macro33(*this); 
    case  34: return macro34(*this); 
    case  35: return macro35(*this); 
    case  36: return macro36(*this); 
    case  37: return macro37(*this); 
    case  38: return macro38(*this); 
    case  39: return macro39(*this); 
    case  40: return macro40(*this); 
    case  41: return macro41(*this); 
    case  42: return macro42(*this); 
    case  43: return macro43(*this); 
    case  44: return macro44(*this); 
    case  45: return macro45(*this); 
    case  46: return macro46(*this); 
    case  47: return macro47(*this); 
    case  48: return macro48(*this); 
    case  49: return macro49(*this); 
    case  50: return macro50(*this); 
    case  51: return macro51(*this); 
    case  52: return macro52(*this); 
    case  53: return macro53(*this); 
    case  54: return macro54(*this); 
    case  55: return macro55(*this); 
    case  56: return macro56(*this); 
    case  57: return macro57(*this); 
    case  58: return macro58(*this); 
    case  59: return macro59(*this); 
    case  60: return macro60(*this); 
    case  61: return macro61(*this); 
    case  62: return macro62(*this); 
    case  63: return macro63(*this); 
    case  64: return macro64(*this); 
    case  65: return macro65(*this); 
    case  66: return macro66(*this); 
    case  67: return macro67(*this); 
    case  68: return macro68(*this); 
    case  69: return macro69(*this); 
    case  70: return macro70(*this); 
    case  71: return macro71(*this); 
    case  72: return macro72(*this); 
    case  73: return macro73(*this); 
    case  74: return macro74(*this); 
    case  75: return macro75(*this); 
    case  76: return macro76(*this); 
    case  77: return macro77(*this); 
    case  78: return macro78(*this); 
    case  79: return macro79(*this); 
    case  80: return macro80(*this); 
    /*
    case  81: return macro81(*this); 
    case  82: return macro82(*this); 
    case  83: return macro83(*this); 
    case  84: return macro84(*this); 
    case  85: return macro85(*this); 
    case  86: return macro86(*this); 
    case  87: return macro87(*this); 
    case  88: return macro88(*this); 
    case  89: return macro89(*this); 
    case  90: return macro90(*this); 
    case  91: return macro91(*this); 
    case  92: return macro92(*this); 
    case  93: return macro93(*this); 
    case  94: return macro94(*this); 
    case  95: return macro95(*this); 
    case  96: return macro96(*this); 
    case  97: return macro97(*this); 
    case  98: return macro98(*this); 
    case  99: return macro99(*this);
    case 100: return macro100(*this); 
    case 101: return macro101(*this); 
    case 102: return macro102(*this); 
    case 103: return macro103(*this); 
    case 104: return macro104(*this); 
    case 105: return macro105(*this); 
    case 106: return macro106(*this); 
    case 107: return macro107(*this); 
    case 108: return macro108(*this); 
    case 109: return macro109(*this); 
    case 110: return macro110(*this); 
    case 111: return macro111(*this); 
    case 112: return macro112(*this); 
    case 113: return macro113(*this); 
    case 114: return macro114(*this); 
    case 115: return macro115(*this); 
    case 116: return macro116(*this); 
    case 117: return macro117(*this); 
    case 118: return macro118(*this); 
    case 119: return macro119(*this); 
    case 120: return macro120(*this); 
    case 121: return macro121(*this); 
    case 122: return macro122(*this); 
    case 123: return macro123(*this); 
    case 124: return macro124(*this); 
    case 125: return macro125(*this); 
    case 126: return macro126(*this); 
    case 127: return macro127(*this); 
    case 128: return macro128(*this); 
    case 129: return macro129(*this); 
    case 130: return macro130(*this); 
    case 131: return macro131(*this); 
    case 132: return macro132(*this); 
    case 133: return macro133(*this); 
    case 134: return macro134(*this); 
    case 135: return macro135(*this); 
    case 136: return macro136(*this); 
    case 137: return macro137(*this); 
    case 138: return macro138(*this); 
    case 139: return macro139(*this); 
    case 140: return macro140(*this); 
    case 141: return macro141(*this); 
    case 142: return macro142(*this); 
    case 143: return macro143(*this); 
    case 144: return macro144(*this); 
    case 145: return macro145(*this); 
    case 146: return macro146(*this); 
    case 147: return macro147(*this); 
    case 148: return macro148(*this); 
    case 149: return macro149(*this); 
    case 150: return macro150(*this); 
    case 151: return macro151(*this); 
    case 152: return macro152(*this); 
    case 153: return macro153(*this); 
    case 154: return macro154(*this); 
    case 155: return macro155(*this); 
    case 156: return macro156(*this); 
    case 157: return macro157(*this); 
    case 158: return macro158(*this); 
    case 159: return macro159(*this); 
    case 160: return macro160(*this); 
    case 161: return macro161(*this); 
    case 162: return macro162(*this); 
    case 163: return macro163(*this); 
    case 164: return macro164(*this); 
    case 165: return macro165(*this); 
    case 166: return macro166(*this); 
    case 167: return macro167(*this); 
    case 168: return macro168(*this); 
    case 169: return macro169(*this); 
    case 170: return macro170(*this); 
    case 171: return macro171(*this); 
    case 172: return macro172(*this); 
    case 173: return macro173(*this); 
    case 174: return macro174(*this); 
    case 175: return macro175(*this); 
    case 176: return macro176(*this); 
    case 177: return macro177(*this); 
    case 178: return macro178(*this); 
    case 179: return macro179(*this); 
    case 180: return macro180(*this); 
    case 181: return macro181(*this); 
    case 182: return macro182(*this); 
    case 183: return macro183(*this); 
    case 184: return macro184(*this); 
    case 185: return macro185(*this); 
    case 186: return macro186(*this); 
    case 187: return macro187(*this); 
    case 188: return macro188(*this); 
    case 189: return macro189(*this); 
    case 190: return macro190(*this); 
    case 191: return macro191(*this); 
    case 192: return macro192(*this); 
    case 193: return macro193(*this); 
    case 194: return macro194(*this); 
    case 195: return macro195(*this); 
    case 196: return macro196(*this); 
    case 197: return macro197(*this); 
    case 198: return macro198(*this); 
    case 199: return macro199(*this);
    case 200: return macro200(*this);
    */
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







