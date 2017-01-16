
#ifndef incl_ArrayFileIO_h
#define incl_ArrayFileIO_h


#include <iostream>
#include <fstream>



template<typename Type> bool readArrayFromFile(char *fileName, Type **x, int *n)
{
  ifstream inFile;

  int i;

  Type dum;

  inFile.open(fileName);

  if (!inFile) return false;

  *n = -1; while (inFile) { inFile >> dum; (*n)++; }

  *x = new Type [*n];

  inFile.clear();
  inFile.seekg(0);

  i = 0; while (i<*n) inFile >> (*x)[i++]; 

  inFile.close();

  return true;
}







template<typename Type> bool writeArrayToFile(char *fileName, Type *x, int n)
{
  ofstream outFile;

  int i;

  outFile.open(fileName);

  if (!outFile) return false;

  i = 0; while (i<n) outFile << x[i++] << "\n"; 

  outFile.close();

  return true;
}





#endif




