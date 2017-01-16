
#include <fstream>
#include <iostream>

#include "MyString.h"
#include "FunctionsProgram.h"
#include "DomainTree.h"


extern DomainTree domain;


bool prgNextDataKey(std::ifstream &Ifile, MyString &keyStrg)
{
  while (Ifile)
    {
       if (keyStrg.begins("BEGIN ")) prgError(1,"prgNextDataKey","'END ...' not found!");
      
       if (domain.isCorrectEndString(keyStrg)) break;
       
       if (keyStrg.begins("END ")) prgError(1,"prgNextDataKey","wrong expression 'END ... ..'!");
      
       if (keyStrg.length() > 0) break;
       
       keyStrg.read(Ifile).stripToMin();
    } 

  if (keyStrg.begins("END ")) return false;
	
  if (!Ifile) return false; 

  return true;
}


