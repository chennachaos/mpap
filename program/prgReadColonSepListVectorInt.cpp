
#include <iostream>

#include "MyString.h"
#include "List.h"
#include "MathVector.h"


// read colon separated list of integer arrays


bool prgReadColonSepListVectorInt(std::ifstream &Ifile, MyString &line, List< Vector<int> > &list)
{
  int i, val, nw = 0, c = 0;

  MyString tmpl, *word;

  list.free();

  while (1)
  {
    list.add(new Vector<int>);

    while (1)
    {
      if (nw == 0)
      {
        line.getNextLine(Ifile);

        nw = line.split(&word);

        i = 0;
      }
      while (i < nw && word[i].toInt(&val,false)) { list[c].append(val); i++; }
    
      if (i < nw) break;

      for (i=0; i<nw; i++) word[i].free(); delete [] word; nw = 0;
    }
    if (word[i] != ":") break; 
    else 
    {
      i++;
      c++;
    }
  }

  //for (i=0; i<list.n; i++) std::cout << list[i] << "\n";

  //exit(0);

  if (i != 0) { list.free(); return false; }

  return true;
}



