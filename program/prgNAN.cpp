

#include "MyString.h"


bool prgNAN(double x)
{
  char tmp[20];

  int i;

  double y;

  sprintf(tmp,"%12.5g  ",x);

  //std::cout << tmp << "= tmp\n\n";

  while (tmp[0] == ' ') { for (i=0; i<strlen(tmp); i++) tmp[i] = tmp[i+1]; }

  while (tmp[strlen(tmp)-1] == ' ') { tmp[strlen(tmp)-1] = '\0'; }

  if (!scanDbl(tmp,&y,false)) return true;
  
  return false;
}







bool prgNAN(double *x, int n)
{
  for (int i=0; i<n; i++) if (prgNAN(x[i])) return true;

  return false;
}









