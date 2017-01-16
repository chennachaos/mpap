
#include <iostream>
#include <fstream>

#include "MyString.h"
#include "FunctionsProgram.h"


using namespace std;


void prgResolveIncludeFile(std::ofstream &tmpFile, char *IfileName)
{
  //cout << IfileName << "\n\n";

  char c;

  int l;
  
  bool comment = false;
  
  MyString line;
  
  ifstream inFile;

  inFile.open(IfileName);

  if (!inFile) prgError(1,"prgResolveIncludeFile","include file not found!'");  
  
  while (inFile) 
  {
    inFile.get(c);

    if (charInCharArray(c,COMMENT_MARK)) comment = true;
    else if (c == '\n')                  comment = false;
    
    if (c == '#' && !comment)
    {
      line.read(inFile);
      line.stripToMin();
      
      if (line.begins("include<"))
      {
	l = 8;
	while (l<line.length() && (line.asCharArray())[l] != '>') l++;
	if (l<line.length())
	{
          line.trunc(l);
	  //cout << &((line.asCharArray())[8]) << " ---> prgResolve....\n\n";
          prgResolveIncludeFile(tmpFile,&(line.asCharArray()[8]));
	}
	else
          prgError(2,"prgResolveIncludeFile","bad statement '#include<File> !'");
      }
      else
	prgError(3,"prgResolveIncludeFile","bad statement '#include<File> !'");

    }
    else
	    
      tmpFile.put(c);
  }

  //cout << "close " << IfileName << "\n\n";
  inFile.close();
  
  return;  
}

