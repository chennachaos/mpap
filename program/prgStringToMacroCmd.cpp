
#include <iostream>

#include "MyString.h"
#include "MacroCommand.h"
#include "MacroList.h"
#include "MacroDialogData.h"
#include "DomainTree.h"
#include "MathBasic.h"
#include "RunControl.h"
#include "FunctionsProgram.h"


extern MacroList  macro;
extern DomainTree domain;
extern RunControl runCtrl;


using namespace std;


// this function translates a string into a macro command and takes 
// into account the default values


bool prgStringToMacroCmd(MacroCommand &macCmd, MyString &macStrg) 
{
   char fct[] = "prgStringToMacroCmd";

   macCmd.free();
	
   //if (runCtrl.mode == NOPROJ) { macCmd.ii = -1; return false; }
	
   char **word = NULL;

   double dp;
   
   int i, j, mi, iw, pc = 0, nw = macStrg.split(&word), cDflt = 0;
   
   MacroDialogData *db;
   
   if (nw < 1) goto rubbish;
   
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// macro type and name
   
   iw = 0; mi = 0; while (mi < macro.n && macro[mi].name != word[iw]) mi++;

   if (mi == macro.n) 
   {  
     if (nw < 2) goto rubbish;
     else 
     {  
       iw = 1; mi = 0; while (mi < macro.n && macro[mi].name != word[iw]) mi++;
       if (mi < macro.n && macro[mi].type != word[iw-1]) goto rubbish; 
     } 
   }

   if (mi == macro.n) goto rubbish;

   macCmd.ii = mi; iw++;

   //cout << macro[macCmd.ii].name << "\n\n";
   
   db = &(macro[mi].db);
   
   // check mode

   if (!macro[mi].sensitivity[runCtrl.mode]) goto wrongMode;
   
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// strg
  
   if (db->strgList.n > 0)
   {  
     if (iw >= nw)  macCmd.strg = db->dfltStrg;
     else
     {
        i = 0; while (i < db->strgList.n) {  if (db->strgList[i] == word[iw]) break;  i++;  }
     
        if (i == db->strgList.n) 
        { 
           if (strlen(word[iw])>0) goto rubbish;
           macCmd.strg = db->dfltStrg; 
        }
        else macCmd.strg = db->strgList[i];
     }
   }
   
   else if (!!db->strgTxtFLabl)
   {
     if (iw >= nw) macCmd.strg = db->dfltStrg;
     else          macCmd.strg = word[iw];
   }
 
   iw++;

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// domain selection

   if (db->selDm)
   {
     if (domain.ndom < 1)
     {
       macCmd.p[pc].x = 0.; pc++; iw++;
       macCmd.p[pc].x = 0.; pc++; iw++;
     }
     else
     {
       i = 0; 
       if (iw < nw && strlen(word[iw]) > 0) 
         if (!macCmd.p[pc].interprete(word[iw],NON_NEG_INT)) goto rubbish;
         else i = roundToInt(macCmd.p[pc].x);

       if (!macCmd.p[pc].fctName)
       {
         if (domain.nDomainOfType(i) < 1) { i = 0; while (domain.nDomainOfType(i) < 1) i++; }
         macCmd.p[pc].x = (double) i;
       }
       pc++; iw++;

       j = 1; 
       if (iw < nw && strlen(word[iw]) > 0)
         if (!macCmd.p[pc].interprete(word[iw],NON_NEG_INT)) goto rubbish;
         else j = roundToInt(macCmd.p[pc].x);
       if (!macCmd.p[pc].fctName)
       {
         if (domain.nDomainOfType(i) < j) j = domain.nDomainOfType(i);
         macCmd.p[pc].x = (double) j;
       }
       pc++; iw++;
     }
   }
   
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// p

   for (i=0; i<db->inputType.dim(); i++)
   {
     switch (db->inputType[i])
     {
       case LIST: if (iw >= nw) { macCmd.p[pc++].x = db->dflt[cDflt]; break; }
	          
                  if (!macCmd.p[pc++].interprete(word[iw++],NON_NEG_INT))
		  {
		    if (strlen(word[iw-1]) > 0) goto rubbish;
		    macCmd.p[pc-1].x = db->dflt[cDflt];
		  }

		  break;
	   
       case RBOX: if (iw >= nw) { macCmd.p[pc++].x = db->dflt[cDflt]; break; }
	          
		  if (!macCmd.p[pc++].interprete(word[iw++],NON_NEG_INT))
		  {
		    if (strlen(word[iw-1]) > 0) goto rubbish;
		    macCmd.p[pc-1].x = db->dflt[cDflt];
		  }
		  
		  break;

       case TBTN: if (iw >= nw) { macCmd.p[pc++].x = db->dflt[cDflt]; break; }
					
		  if (!macCmd.p[pc++].interprete(word[iw++],NON_NEG_INT))
		  {
		    if (strlen(word[iw-1]) > 0) goto rubbish;
		    macCmd.p[pc-1].x = db->dflt[cDflt];
		  }
		  
		  break;
		  
       case TXTF: if (iw >= nw) { macCmd.p[pc++].x = db->dflt[cDflt]; break; }
                
		  if (!macCmd.p[pc++].interprete(word[iw++],ALL_REAL))
		  {
		    if (strlen(word[iw-1]) > 0) goto rubbish;
	            macCmd.p[pc-1].x = db->dflt[cDflt]; 
		  }
		  
                  break;
     }

     if (db->inputType[i] != LABL) cDflt++;
   }   

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
   
   //quit:

   i = 0; if (word != NULL) while (word[i] != NULL) delete [] word[i++]; delete [] word;

   if (macCmd.ii < 0) prgError(1,"prgStringToMacroCmd","macCmd.ii < 0!");
     
   return true;

   rubbish:

     macCmd.ii = -1; 
     
     i = 0; if (word != NULL) while (word[i] != NULL) delete [] word[i++]; delete [] word;

     std::cout << "\n         '" << macStrg << "' is RUBBISH!\n\n";

     return false;

   wrongMode:

     macCmd.ii = -1; 

     i = 0; if (word != NULL) while (word[i] != NULL) delete [] word[i++]; delete [] word;

     std::cout << "\n         '" << macStrg << "' is not available in this mode!\n\n";
     
     return false;
}


