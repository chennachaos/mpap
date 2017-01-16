
#include <cmath>

#include "MacroCommand.h"
#include "MacroList.h"
#include "MacroDialogData.h"


extern MacroList  macro;


// this function translates the macro command to a string without 
// any consideration of default values


void prgMacroCmdToSimpleString(MyString &macCmdStr, MacroCommand &macCmd)
{
  char tmp[50];

  int i, n;

  MacroDialogData &db = macro[macCmd.ii].db;
	
  macCmdStr.free();

  macCmdStr.append(macro[macCmd.ii].type).append(',').append(macro[macCmd.ii].name);

  if (db.strgList.n > 0 || !!db.strgTxtFLabl)  macCmdStr.append(',').append(macCmd.strg);
  
  else if (db.inputType.dim() > 0 || db.selDm) macCmdStr.append(',');

  n = db.inputType.dim() - db.labl.n;  if (db.selDm) n = n + 2;
  
  for (i=0; i<n; i++) 
  {
    if (!macCmd.p[i].fctName) 
    {  
      sprintf(tmp,"%g",macCmd.p[i].x);
      macCmdStr.append(',').append(tmp);
    }
    else
    {
      macCmdStr.append(",$").append(macCmd.p[i].fctName);
      if (!!macCmd.p[i].fctPar) macCmdStr.append('(').append(macCmd.p[i].fctPar).append(')');
    }
  }
  
  return;  
}

