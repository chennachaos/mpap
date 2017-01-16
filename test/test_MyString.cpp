
#include <iostream>

#include "MyString.h"


bool debug = true;

int main()
{
   cout << "\n  Hello, world.\n\n";

   MyString tmp;

   cout << tmp << "|\n\n";
   
   while (1)
   {
      cout << " input some stuff : ";
      
      tmp.read(cin);
      
      cout << "\n you have typed: <" << tmp << ">\n\n";

      if (tmp.stripToMin().begins("quit")) break;
   }
   
   return 0;
};




