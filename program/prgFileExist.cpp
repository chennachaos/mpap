
#include <iostream>
#include <fstream>

#include "MyString.h"
#include "Definitions.h"
#include "Debug.h"




bool prgFileExist(char *directory, char *file)
{
  if (file == NULL || strlen(file) < 1) return false;
  
  int l2 = strlen(file), l1 = 0;

  if (directory != NULL) l1 = strlen(directory);
  
  char *pathAndFile = new char[l1 + l2 + 2];
  
  if (directory == NULL)  sprintf(pathAndFile,file);
    else sprintf(pathAndFile,"%s%c%s",directory,SLASH,file);
  
  if (debug) std::cout << " " << pathAndFile << "\n\n";

  std::ifstream stream;

  stream.open(pathAndFile);

  delete [] pathAndFile;
  
  if (!stream) return false;

  stream.close();

  return true;	
}





bool prgFileExist(char *pathAndFile)
{
  if (pathAndFile == NULL) return false;

  if (strlen(pathAndFile) < 1) return false;
	
  std::ifstream stream;

  stream.open(pathAndFile);

  if (!stream) return false;

  stream.close();

  return true;	
}







bool prgFileExist(MyString &directory, MyString &file)
{
  return prgFileExist(directory.asCharArray(),file.asCharArray());
}



bool prgFileExist(MyString &pathAndFile)
{
  return prgFileExist(pathAndFile.asCharArray());
}


