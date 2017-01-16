
#include <iostream>

#include "MacroList.h"
#include "Definitions.h"
#include "FunctionsProgram.h"
#include "Debug.h"


using namespace std;


//  MacroList constructor

MacroList::MacroList(void) 
{ 
  ntype = 0;

  jp = NULL;

  if (debug) cout << " MacroList constructor \n\n";
}






// MacroList destructor

MacroList::~MacroList() 
{ 
  delete [] jp;

  if (debug) cout << " MacroList destructor \n\n";
}
 






//  check Macro Names

int MacroList::checkNames(void)
{
  int i, j;

  MacroList &macro = *this;

  for (i=0; i<macro.n; i++) 
  { 
    if (macro[i].name.length()<1 || macro[i].type.length()<1) return macro[i].id;
    
    for (j=0; j<macro.n; j++)   if (macro[i].type == macro[j].name) return macro[j].id;
    
    for (j=i+1; j<macro.n; j++) if (macro[i].name == macro[j].name) return macro[j].id;
  }

  return -1;
}







//  sort Macros according to macro type

void MacroList::sort(void)
{    
  int i = 0, j;

  MacroList &macro = *this;
  
  if (n == 0) return;

  ntype = 0;
 
  do
  {  j = i;
     do   if (macro[i].type == macro[++j].type)  {  i++; if (i != j) move(j,i); }
     while (j < n - 1);
     if ( ++i<n && macro[i].type != macro[i-1].type)  ntype++; 
  }  while (i < n - 1);
  ntype++; 

  // calculate profile jp

  jp = new int[ntype];

  jp[ntype-1] = n - 1;
  
  j = 0;
  for (i=0;i<ntype-1;i++) 
    {  while (macro[j].type == macro[j+1].type) { j++; }  jp[i] = j++;  }
 
  if (debug) for (i=0;i<ntype;i++) { cout << " " << jp[i]; cout << "\n\n"; }
  
  return;
}









