
#include "Macro.h"
#include "List.h"
#include "Counter.h"
#include "FunctionsProgram.h"


extern ListInfinite<Counter> counter;


int macro42(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "coun";
    macro.type = "ctrl";
    macro.what = "increment and print counter";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;
    macro.sensitivity[PRE]   = true;

    macro.db.stringTextField("name :","",15);
  
    macro.db.addTextField("",1.); 
   
    macro.db.addList("set","*incr","mult"); 

    macro.db.addToggleButton("suppress output",false);

    macro.db.addToggleButton("write to Tfile",false);
    
    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int i = roundToInt(macro.p[0]),
      j = roundToInt(macro.p[1]), k;

  MyString name;

  char tmp[30];

  name = macro.strg;

  k = 0; while (k<counter.n && counter[k].name != name) k++;
  if (k==counter.n && j>1)
    { prgWarning(1,"macro42","'ctrl,coun' ignored! invalid counter name!"); return 0; }

  switch (j)
  {
    case  1: if (name.length() < 1 || !name) 
	       { prgWarning(1,"macro42","'ctrl,coun' ignored! invalid counter name!"); return 0; }
	     
	     if (k==counter.n) counter[k].name = name;
	     
	     counter[k].val = i;  break;  // set to value
	    
    case  2: counter[k].val += i; break;  // increment counter

    case  3: counter[k].val *= i; break;  // multiply counter
 
    default: prgWarning(2,"macro42","'ctrl,coun' ignored! unknown operation id!"); return 0;
  }

  // print on screen
  
  if (roundToInt(macro.p[2]) == 0)
    COUT << "counter(" << counter[k].name << ") = " << counter[k].val << "\n\n";

  // write to Tfile

  if (roundToInt(macro.p[3]) == 1)
  {
    sprintf(tmp,"%d",counter[k].val);

    prgWriteToTFile(tmp);
  }
 
//--------------------------------------------------------------------------------------------------
  return 0;  
}

