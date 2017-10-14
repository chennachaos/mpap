
#ifndef incl_MyString_h
#define incl_MyString_h

#include <iostream>
#include <fstream>
#include <string.h>
#include <cstring>

#include "List.h"
#include "MathVector.h"


class MyString: public ListItem
{
  public:

    MyString(void);
    MyString(char*);
    virtual ~MyString();
    
    char &operator[](int);

    bool operator==(MyString &);
    bool operator!=(MyString &);

    bool operator==(char*);
    bool operator!=(char*);
    
    bool operator!(void);

    bool begins(char*);
    bool begins(MyString &);
   
    bool contains(char*);
    bool contains(MyString &);
    
    int containsAtPos(char*);
    int containsAtPos(MyString &);
    
    MyString &operator=(MyString &);
    MyString &operator=(char *);

    MyString &append(MyString &);
    MyString &append(char*);
    MyString &append(char);

    MyString &insert(int, char*);
    MyString &insert(int, MyString &);
    
    MyString &trunc(int);

    MyString &free(void);
    
    char* asCharArray(void);
    
    int length(void);
    
    MyString &read(std::istream &);

    MyString &getNextLine(std::istream &);
    
    MyString &input(void);
    MyString &inputKeepIfReturn(void);

    void print(void);

    int which(char**);
    int which(MyString*);

    bool FileOK(void);
    bool DirectoryOK(void);
    
    MyString &strip(void);
    MyString &stripToMin(void);

    int countWords(void);
  
    bool copyAfter(char, MyString &);
    
    int split(char***);
    int split(MyString**);
    
    bool toInt(int *, bool errMsg = true);
    bool toDbl(double *, bool errMsg = true);

    bool toInt(Vector<int> &, bool errMsg = true);

  private:

    char* ch;
};



MyString *MyStringArray(char**);
std::ostream &operator<<(std::ostream  &, MyString &);
std::ofstream &operator<<(std::ofstream &, MyString &);




inline bool MyString::operator!(void)
{
  if (ch == NULL) return true; 

  if (length() == 0) return true;

  return false;
}




inline MyString &MyString::operator=(MyString &strg2)
{
  free();   append(strg2);   return *this;
}



inline MyString &MyString::operator=(char *strg2)
{
  free();   append(strg2);   return *this;
}



inline bool charInCharArray(char c, char *s)
{
  int i = 0, l = strlen(s); 
  
  while(i < l && c != s[i]) i++;

  if (i < l) return true;

  return false;
}





inline bool scanInt(char *s, int *i, bool errMsg = true)
{
  int j = 0;
  if (s == NULL) goto error;
  
  if (strlen(s) == 0) goto error;
  
  while (charInCharArray(s[j],"1234567890+-")) j++;  
  if (s[j]!='\0') goto error;
  
  if (!sscanf(s,"%d",i)) goto error;
  
  return true;
  
  error:
    if (errMsg) std::cout << " ERROR in 'scanInt'!";   
    return false;
}  






inline bool scanDbl(char *s, double *x, bool errMsg = true)
{
  int j = 0; bool dot = false, E = false;	
  if (s == NULL) goto error;

  if (strlen(s) == 0) goto error;
  
  while (charInCharArray(s[j],"1234567890.-+eE")) 
  { 
    if (s[j] == '.')                 { if (!dot) dot = true; else goto error; }
    else if (s[j]=='e' || s[j]=='E') { if (!E)     E = true; else goto error; }
    j++;
  }
  if (s[j]!='\0') goto error;
  
  if (!sscanf(s,"%lf",x)) goto error;
  
  return true;
  
  error:
    if (errMsg) std::cout << " ERROR in 'scanDbl'!";   
    return false;
}  





inline int which(char *word, char **list)
{
  if (word == NULL || list == NULL) return -1;

  int i = 0;   while (list[i] != NULL && strcmp(word,list[i]) != 0) i++;

  if (list[i] == NULL) return -1;
  
  return i;
}






#define ALL_WORD_SEPARATORS  " ,;!|*?\t"
#define WEAK_WORD_SEPARATORS " "

#define MAX_READ_LINE_LENGTH 100

#define COMMENT_MARK "!%"

#define APPEND_NEXT_LINE_MARK '&'

#endif



