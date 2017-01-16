
#include "FunctionsProgram.h"
#include "List.h"
#include "MyStringList.h"
#include "Debug.h"



MyStringList::MyStringList(void)
{
  if (debug) std::cout << " constructor MyStringList\n\n";
}





MyStringList::~MyStringList()
{
  if (debug) std::cout << " destructor MyStringList\n\n";
}





void MyStringList::addNew(char* s0, 
                          char* s1 , char* s2 , char* s3 , char* s4 , 
                          char* s5 , char* s6 , char* s7 , char* s8 , 
                          char* s9 , char* s10, char* s11, char* s12, 
                          char* s13, char* s14, char* s15, char* s16, 
                          char* s17, char* s18, char* s19, char* s20, 
                          char* s21, char* s22, char* s23, char* s24, 
                          char* s25, char* s26, char* s27, char* s28, 
                          char* s29, char* s30, char* s31, char* s32, 
                          char* s33, char* s34, char* s35, char* s36, 
                          char* s37, char* s38, char* s39, char* s40 )
{
  if (s0  != NULL) add(new MyString(s0 ));
  if (s1  != NULL) add(new MyString(s1 ));
  if (s2  != NULL) add(new MyString(s2 ));
  if (s3  != NULL) add(new MyString(s3 ));
  if (s4  != NULL) add(new MyString(s4 ));
  if (s5  != NULL) add(new MyString(s5 ));
  if (s6  != NULL) add(new MyString(s6 ));
  if (s7  != NULL) add(new MyString(s7 ));
  if (s8  != NULL) add(new MyString(s8 ));
  if (s9  != NULL) add(new MyString(s9 ));
  if (s10 != NULL) add(new MyString(s10));
  if (s11 != NULL) add(new MyString(s11));
  if (s12 != NULL) add(new MyString(s12));
  if (s13 != NULL) add(new MyString(s13));
  if (s14 != NULL) add(new MyString(s14));
  if (s15 != NULL) add(new MyString(s15));
  if (s16 != NULL) add(new MyString(s16));
  if (s17 != NULL) add(new MyString(s17));
  if (s18 != NULL) add(new MyString(s18));
  if (s19 != NULL) add(new MyString(s19));
  if (s20 != NULL) add(new MyString(s20));
  if (s21 != NULL) add(new MyString(s21));
  if (s22 != NULL) add(new MyString(s22));
  if (s23 != NULL) add(new MyString(s23));
  if (s24 != NULL) add(new MyString(s24));
  if (s25 != NULL) add(new MyString(s25));
  if (s26 != NULL) add(new MyString(s26));
  if (s27 != NULL) add(new MyString(s27));
  if (s28 != NULL) add(new MyString(s28));
  if (s29 != NULL) add(new MyString(s29));
  if (s30 != NULL) add(new MyString(s30));
  if (s31 != NULL) add(new MyString(s31));
  if (s32 != NULL) add(new MyString(s32));
  if (s33 != NULL) add(new MyString(s33));
  if (s34 != NULL) add(new MyString(s34));
  if (s35 != NULL) add(new MyString(s35));
  if (s36 != NULL) add(new MyString(s36));
  if (s37 != NULL) add(new MyString(s37));
  if (s38 != NULL) add(new MyString(s38));
  if (s39 != NULL) add(new MyString(s39));
  if (s40 != NULL) add(new MyString(s40));
      
  return;
}





void MyStringList::print(void)
{
  for (int i=0; i<n; i++)  std::cout << " " << (*this)[i] << "\n"; std::cout << "\n";  return;
}





char **MyStringList::genCharArrays(void)
{
  char **tmp = new char* [n+1];

  for (int i=0; i<n; i++)  
    {  tmp[i] = new char[(*this)[i].length()];  sprintf(tmp[i],(*this)[i].asCharArray());  }

  tmp[n] = NULL;

  return tmp;
}







int MyStringList::which(char *strg)
{
  int i = 0;
  
  while (i<this->n && (*this)[i] != strg) i++;

  if (i == this->n) return -1;
  
  return i;
}




int MyStringList::which(MyString &strg)
{
  int i = 0;
  
  while (i<this->n && (*this)[i] != strg) i++;

  if (i == this->n) return -1;
  
  return i;
}






int MyStringList::whichBegins(char *strg)
{
  int i = 0;
  
  while (i<this->n && strncmp(strg,(*this)[i].asCharArray(),strlen(strg))!=0)  i++;

  if (i == this->n) return -1;
  
  return i;
}




int MyStringList::whichBegins(MyString &strg)
{
  int i = 0;
 
  while (i<this->n && !strg.begins((*this)[i])) i++;

  if (i == this->n) return -1;
  
  return i;
}









// overload printing to screen

std::ostream &operator<<(std::ostream &os,MyStringList &list)
{
  os << "{";
  for (int i=0; i<list.n-1; i++) os << list[i] << ","; 
  if (list.n > 0) os << list[list.n-1] << "}"; else os << " }";

  return os;
}


