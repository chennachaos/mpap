
#include <iostream>


#include "PropertyItem.h"
#include "FunctionsProgram.h"
#include "Debug.h"


using namespace std;



PropertyItem::PropertyItem(void)
{
  return;
}



PropertyItem::PropertyItem(int propType)
{
  type = propType;

  if (debug) cout << " PropertyItem constructor\n\n";
  
  return;
}






PropertyItem::~PropertyItem()
{
  if (debug) cout << " PropertyItem destructor\n\n";
  
  return;
}







void PropertyItem::readInputData(std::ifstream &Ifile, MyString &line, char *errMsg)
{
  int nw, i, j;

  char fct[] = "PropertyItem::readInputData";

  myVector<double> dataTmp;

  MyString *word;
  
  line.getNextLine(Ifile);
 
  nw = line.split(&word);
 
  if (nw < 2)              prgError(1,fct,errMsg);

  if (!word[0].toInt(&id)) prgError(2,fct,errMsg);
  
  name = word[1];
  
  for (j=2; j<nw; j++)   if (!word[j].toDbl(dataTmp.append()))  prgError(4,fct,errMsg);
    
  while (1)
  {
    for (i=0; i<nw; i++) word[i].free(); delete [] word;
    
    line.getNextLine(Ifile); 

    nw = line.split(&word);

    j = 0;
    
    while (j<nw && word[j].toDbl(dataTmp.append(),false))  j++;

    if (j<nw) { for (i=0; i<j+1; i++) dataTmp.del(dataTmp.n-1); break; }
  }

  if (dataTmp.n == 0) dataTmp.append(0);

  data = dataTmp;

  for (i=0; i<nw; i++) word[i].free(); delete [] word;
  
  if (j>0 && j<nw) prgError(5,fct,errMsg);

  return;
}






