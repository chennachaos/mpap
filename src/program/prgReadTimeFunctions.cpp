
#include <iostream>
#include <fstream>

#include "MyString.h"
#include "TimeFunction.h"
#include "FunctionsProgram.h"
#include "Debug.h"


extern List<TimeFunction> timeFunction;

using namespace std;



void prgReadTimeFunctions(std::ifstream &Ifile)
{
  MyString line, *word;

  int i, ii, j, j1, id, nw;

  line.getNextLine(Ifile).stripToMin();  

  while (Ifile && line != "END TIME_FUNCTIONS")
  {
    nw = line.split(&word);

    if (nw < 4) prgError(1,"prgReadTimeFunctions","error in data!");

    if (!word[0].toInt(&id)) prgError(2,"prgReadTimeFunctions","error in data!");

    if (id == 0) prgError(2,"prgReadTimeFunctions","0 as id not admissible!");

    i = 0; while (i<timeFunction.n && id!=timeFunction[i].id) i++;

    if (i == timeFunction.n)  { timeFunction.add(new TimeFunction); timeFunction[i].id = id; }

    timeFunction[i].fct.add(new TimeFunctionCore);
    ii = timeFunction[i].fct.n - 1;

    if (!word[1].toDbl(&timeFunction[i].fct[ii].t0)) 
      prgError(3,"prgReadTimeFunctions","error in data!");

    if (!word[2].toDbl(&timeFunction[i].fct[ii].t1)) 
      prgError(4,"prgReadTimeFunctions","error in data!");

    j1 = nw; if (nw>11) j1 = 11;

    for (j=3; j<j1; j++) 
      if (!word[j].toDbl(&timeFunction[i].fct[ii].p[j-3])) 
        prgError(5,"prgReadTimeFunctions","error in data!");

    for (i=0; i<nw; i++) word[i].free(); delete [] word;

    line.getNextLine(Ifile).stripToMin();  
  }

  for (i=0; i<timeFunction.n; i++) timeFunction[i].update();

  if (debug)
    { for (i=0; i<timeFunction.n; i++)   cout << timeFunction[i] << "\n";  cout << "\n"; }

  return;	
}



