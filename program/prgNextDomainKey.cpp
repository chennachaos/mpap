
#include <iostream>
#include <fstream>

#include "DomainTree.h"
#include "MyString.h"
#include "FunctionsProgram.h"
#include "Debug.h"


extern DomainTree domain;


int prgNextDomainKey(std::ifstream &Ifile)
{
   char     tmp1[10], *domKey[] = DOMAIN_KEY;

   MyString tmp2, line;
   
   int      i, nDomType = prgListCharArrayLength(domKey);
   
   line.read(Ifile).stripToMin();
   while (Ifile)
     {
	if (line == "BEGIN RUN_CONTROL") { i = -1; break; }
	
	if (line == "BEGIN TIME_FUNCTIONS") { i = -2; break; }
		
	i = 0;
	while (i < nDomType)
	  {	
             tmp2.free().append("BEGIN ").append(domKey[i]).append(' ');

	     if (line.begins(tmp2))
	       {   
                  sprintf(tmp1,"%i",domain.nDomainOfType(i)+1);
                  tmp2.append(tmp1);
		  if (tmp2 == line) break; 
		  std::cout << "   SERIOUS WARNING! wrong identifier in " << tmp2 << "\n\n";
	       }
	     i++;
	  }	
	     
	if (tmp2 == line) break;
	
        line.read(Ifile).stripToMin();
     } 

   if (!Ifile) return -3; 

   return i;
}


