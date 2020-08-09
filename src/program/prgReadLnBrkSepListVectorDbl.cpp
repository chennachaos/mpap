
#include <iostream>

#include "MyString.h"
#include "List.h"
#include "MathVector.h"


// read colon separated list of integer arrays


bool prgReadLnBrkSepListVectorDbl(std::ifstream &Ifile, MyString &line, List< Vector<double> > &list)
{
  int i, j, nw = 0, n = 0;
  double val;

  MyString tmpl, *word;

  list.free();
  
  while (1)
  {
    list.add(new Vector<double>);

    line.getNextLine(Ifile);

    nw = line.split(&word);

    i = 0;  while (i < nw && word[i].toDbl(&val,false)) { list[n].append(val); i++; }

    for (j=0; j<nw; j++) word[j].free(); delete [] word;
    //cout << list[i] << endl;

    if (i < nw) break;

    n++;
  }
  
  list.del(n);

  //for (i=0; i<list.n; i++) std::cout << list[i] << "\n";

  if (n == 0) { list.free(); return false; }

  return true;
}



