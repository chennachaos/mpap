
#include <iostream>

#include "Global.h"
#include "UnixGlobal.h"

#include "MathProfileList.h"

using namespace std;



int main(void)
{
  int i, j;
	
  ProfileListFortran<double> x;

  
  x(1,2) = 1.1;
  x.print();
  
  cout << "-1----------------------------\n\n";
  x(2,1) = 2.2;
  x.print();
  
  cout << "-2----------------------------\n\n";
  x(1,1) = 3.3;
  x.print();
  
  cout << "-3----------------------------\n\n";
  x(2,4) = 4.4;
  x.print();
  
  cout << "-4----------------------------\n\n";
  x(1,4) = 5.5;
  x.print();
  
  cout << "-5----------------------------\n\n";
  x(4,2) = 6.6;
  x.print();
  
  cout << "-6----------------------------\n\n";
  x(3,3) = 7.7;
  x.print();

  cout << "-7----------------------------\n\n";

  for (i=1; i<6; i++)
  {
    for (j=1; j<6; j++)  if (x.nonZero(i,j)) cout << " x"; else cout << " o";
    cout << "\n";
  }
  cout << "\n";
  
  x.zero();

  x.print();
  

  return 0;

};


