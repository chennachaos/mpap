

#include <iostream>


int prgListCharArrayLength(char **list)
{
   int i = 0;

   while (list[i] != NULL) i++;

   return i;
}

