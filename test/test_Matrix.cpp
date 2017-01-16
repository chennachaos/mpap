
#include <iostream>

#include "Global.h"
#include "UnixGlobal.h"

#include "MathMatrix.h"

using namespace std;


int main()
{
  int i, n, m;

  //debug = true;
 
  cout << "Hello, world!\n";
  
  MatrixSparse<double> A;

  A(1,3) = 7; 
  A(3,4) = 2;
  
//  A.printAll();
  
//  A.setMode(ARRAY);

//  A.printAll();

  A.print(cout);

  A(1,1) = 2;
  A(3,2) = 1;
  A(3,3) = 4;

  A.print(cout);

//  A.printAll();

  //A.show(10,2);
 
  
  return 0;
};


