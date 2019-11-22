
#include "MacroDialogData.h"
#include "MyString.h"
#include "FunctionsProgram.h"


MacroDialogData::MacroDialogData(void) 
{ 
   selDm = false; 

   frmRBoxFlg = false;

   frmBBoxFlg = false;

   return; 
}
 

MacroDialogData::~MacroDialogData() 
{ 
   return; 
}
 


void MacroDialogData::selectDomain(void)
{
   selDm = true;  return;
}




bool MacroDialogData::operator!(void)
{
  if (!!strgTxtFLabl || strgList.n > 0 || list.n > 0 || rBox.n > 0 
                         || tBtn.n > 0 || txtF.n > 0 || labl.n > 0) return false;
  
  return true;
}





void MacroDialogData::stringList(char *s1, char *s2, char *s3, char *s4, char *s5, char *s6, 
    	                         char *s7, char *s8, char *s9, char *s10, char *s11, char *s12)
{
  if (s1 == NULL) return;

  if (strgList.n > 0 || !!strgTxtFLabl) prgError(1,"MacroDialogData",
        "You can use only ONE string list or string text field per macro!");
  
  strgList.addNew(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12);

  int d = 0;  
   
  while (d < strgList.n && strgList[d].asCharArray()[0] != '*') d++;
  
  if (d < strgList.n) 
   { 
     dfltStrg = &(strgList[d].asCharArray()[1]); 
     strgList[d].free().append(dfltStrg); 
   }
  else 
     dfltStrg = strgList[0];

  return;   
} 







void MacroDialogData::stringTextField(char *labl, char *s, int c)
{
  if (labl == NULL) return;

  if (strgList.n > 0 || !!strgTxtFLabl) prgError(2,"MacroDialogData",
         "You can use only ONE string list or string text field per macro!");

  strgTxtFLabl = labl;

  strgTxtFCol = c;
  
  dfltStrg = s;

  if (!dfltStrg) dfltStrg = "";
  
  return;
}








void MacroDialogData::addList(char *s1, char *s2, char *s3, char *s4, char *s5, char *s6, 
    	                      char *s7, char *s8, char *s9, char *s10, char *s11, char *s12)
{
   if (s1 == NULL) return;
   
   inputType.append(LIST);
   
   list.add(new MyStringList);
 
   list[list.n-1].addNew(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12); 

   int i = list.n - 1, d = 0;  
   
   while (d < list[i].n && list[i][d].asCharArray()[0] != '*') d++;
  
   MyString tmp;
   
   if (d < list[i].n) { tmp = list[i][d]; list[i][d].free().append(&(tmp.asCharArray()[1])); }
   else d = 0;
   
   dflt.append((double) (d + 1));
   
   return; 
} 




void MacroDialogData::addList(char **s)
{
   if (s == NULL || s[0] == NULL) return;
   
   inputType.append(LIST);
   
   list.add(new MyStringList);
 
   int i = list.n - 1, d = 0;  

   while (s[list[i].n] != NULL) list[i].add(new MyString(s[list[i].n]));
  
   while (d < list[i].n && list[i][d].asCharArray()[0] != '*') d++;
  
   MyString tmp;
   
   if (d < list[i].n) { tmp = list[i][d]; list[i][d].free().append(&(tmp.asCharArray()[1])); }
   else d = 0;
   
   dflt.append((double) (d + 1));
   
   return; 
} 




void MacroDialogData::addRadioBox(char *s1, char *s2, char *s3, char *s4, char *s5, char *s6, 
    	                          char *s7, char *s8, char *s9, char *s10, char *s11, char *s12)
{
   if (s1 == NULL) return;
	
   inputType.append(RBOX);
   
   rBox.add(new MyStringList);
   
   rBox[rBox.n-1].addNew(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12); 
   
   int i = rBox.n - 1, d = 0;  
   
   while (d < rBox[i].n && rBox[i][d].asCharArray()[0] != '*') d++;
  
   MyString tmp;
   
   if (d < rBox[i].n) { tmp = rBox[i][d]; rBox[i][d].free().append(&(tmp.asCharArray()[1])); }
   
   dflt.append((double) (d + 1));
   
   return; 
} 





void MacroDialogData::addLabel(char *s)
{
   if (s == NULL) return;

   autoNextButtonBox();  // if required generate new button box

   nBttn[nBttn.dim()-1]++;  // increase number of buttons in last button box by one
   
   inputType.append(LABL); // append new button type and set to LABL
   
   labl.addNew(s);

   return; 
} 






void MacroDialogData::addToggleButton(char *s, bool dfltVal)
{
   if (s == NULL) return;

   autoNextButtonBox();

   nBttn[nBttn.dim()-1]++;
   
   inputType.append(TBTN);

   dflt.append((double) dfltVal);

   tBtn.addNew(s);
   
   return; 
} 






void MacroDialogData::addTextField(char *s, double x, int c)
{
   if (s == NULL) return;

   autoNextButtonBox();

   nBttn[nBttn.dim()-1]++;
   
   inputType.append(TXTF);
   
   dflt.append(x);
   
   txtF.addNew(s);

   txtFCol.append(c);
   
   return; 
} 




void MacroDialogData::nextButtonBox(void)
{
   if (nBttn.dim() > 0) 
     if (nBttn[nBttn.dim()-1] > 0) 
     {  nBttn.append(0);
        while (bttnCol.n < nBttn.n-1) bttnCol.append(0);
        bttnCol.append(4);
     }
   return;
}
   



void MacroDialogData::autoNextButtonBox(void)
{
   if (inputType.dim() == 0) 
   {  nBttn.append(0); bttnCol.append(4); 
   }
   else 
   {  int i = inputType[inputType.dim()-1];
      if (i == LIST || i == RBOX) 
	if (!(nBttn.dim() > 0 && nBttn[nBttn.dim()-1] == 0)) 
	{  nBttn.append(0);
           while (bttnCol.n < nBttn.n-1) bttnCol.append(0);
           bttnCol.append(4);
	}
   }
   return;
}
      


void MacroDialogData::setButtonColumnDim(int c1, int c2, int c3, int  c4, int  c5, int  c6,
		                         int c7, int c8, int c9, int c10, int c11, int c12)
{
   if (nBttn.dim() >  0) bttnCol[ 0] =  c1; 
   if (nBttn.dim() >  1) bttnCol[ 1] =  c2; 
   if (nBttn.dim() >  2) bttnCol[ 2] =  c3; 
   if (nBttn.dim() >  3) bttnCol[ 3] =  c4; 
   if (nBttn.dim() >  4) bttnCol[ 4] =  c5; 
   if (nBttn.dim() >  5) bttnCol[ 5] =  c6; 
   if (nBttn.dim() >  6) bttnCol[ 6] =  c7; 
   if (nBttn.dim() >  7) bttnCol[ 7] =  c8; 
   if (nBttn.dim() >  8) bttnCol[ 8] =  c9; 
   if (nBttn.dim() >  9) bttnCol[ 9] = c10; 
   if (nBttn.dim() > 10) bttnCol[10] = c11; 
   if (nBttn.dim() > 11) bttnCol[11] = c12; 

   return;
} 





void MacroDialogData::frameRadioBox(void)
{
   frmRBoxFlg = true; return;
}






void MacroDialogData::frameButtonBox(void)
{
   frmBBoxFlg = true; return;
}




