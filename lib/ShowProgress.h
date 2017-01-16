
#ifndef incl_ShowProgress_h
#define incl_ShowProgress_h


#include "MyString.h"
#include "FunctionsProgram.h"



template<typename Type> class ShowProgress
{ 
  public:

    ShowProgress(void) { }

    ~ShowProgress() { }
    
    bool percentageFlg;
    
    int  prevPerc;

    MyString fmt;

    Type total;

    void init(char*, char*, Type, bool);

    void show(Type);
};



template<typename Type> void ShowProgress<Type>::init(char *strg, char *strg2, Type n, bool flg)
{
  char fct[] = "ShowProgress::init";

  percentageFlg = flg;
  
  prevPerc = -1;

  MyString tmp;

  total = n;

  if (percentageFlg) fmt = "%3d"; else fmt = strg2;

  char *ch = fmt.asCharArray();

  int i = 0, l = strlen(ch), d;

  while (i < l) if (ch[i] == '%') break; else i++;

  if (ch[i] != '%' || i == l) prgError(1,fct,"invalid format string!");

  while (!charInCharArray(ch[l-1],"1234567890")) { ch[l-1] = ch[l]; l--; }

  if (!scanInt(ch+i+1,&d,false)) prgError(2,fct,"invalid format string!");

  if (percentageFlg) fmt = "%3d"; else fmt = strg2;

  //cout << d << " = d\n";

  if (percentageFlg)
  {
    tmp.append(strg).append(" ").append(fmt).append("% / ").append(fmt).append("%");

    for (i=0; i<8; i++) tmp.append("\b");

    printf(tmp.asCharArray(),0,100);
  }
  else
  {
    tmp.append(strg).append(" ").append(fmt).append(" / ").append(fmt);

    for (i=0; i<d+3; i++) tmp.append("\b");

    printf(tmp.asCharArray(),0,total);
  }

  for (i=0; i<d; i++) fmt.insert(0,"\b");

  //cout << "|" << tmp << "|\n";
  //cout << "|" << fmt << "|\n";

  return;
}



template<typename Type> void ShowProgress<Type>::show(Type curr)
{
  int perc;

  if (percentageFlg)
  {
    perc = int(curr/total*(Type)100);
    if (perc > prevPerc)
    {
      printf(fmt.asCharArray(),perc);
      prevPerc = perc;
    }
  }
  else
  {
    printf(fmt.asCharArray(),curr);
  }

  if (curr >= total) cout << "\n";

  return;
}











#endif




