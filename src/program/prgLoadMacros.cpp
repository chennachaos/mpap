
#include <iostream>

#include "FunctionsProgram.h"
#include "Definitions.h"
#include "MacroList.h"
#include "Debug.h"


extern MacroList macro;

using namespace std;


void prgLoadMacros(void)

{
  int i, d = 1;

  
  // load all MAX_MACRO macros
	
  for (i=1;i<=MAX_MACRO;i++)  macro.add(new Macro(i));

  if (debug) cout << " macros have been generated\n\n";
  

  // initialise all macros
  
  for (i=1;i<=MAX_MACRO;i++)  { if (debug) cout << i << "\n"; macro[i-1].exec(); }

  if (debug) cout << " macros have been initialised\n\n";

  
  // remove undefined macros
  
  for (i=1;i<=MAX_MACRO;i++)
    if (!macro[i-d])  { macro.del(&macro[i-d]); d++; }

  if (debug) cout << " superfluous macros have been deleted\n\n";

 
  // sort macros according to type
  
  macro.sort();
  
  if (debug) cout << " macros have been sorted\n\n";


  // check names and types

  i = macro.checkNames(); 
  if (i > -1) { cout << " ERROR! type/name of macro " << i << " not admissible!\n\n"; exit(0); }
  
  if (debug) cout << " macro names and types have been checked\n\n";

  
  if (debug) for (i=0;i<macro.n;i++) 
    cout << " macro " << macro[i].id << " (" << macro[i].type << "," << macro[i].name << ")\n\n";
  
  if (debug) cout << " " << macro.n << " macros loaded.\n\n";

  return;
}

