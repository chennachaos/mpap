
#include <iostream>
#include <string>

#include "MyString.h"
#include "Debug.h"


using namespace std;


MyString::MyString(void)
{
  ch = NULL;

  if (debug) std::cout << " MyString constructor\n\n";
  
  //std::cout << ++nObj << "\n\n";
  
  return;
}



MyString::MyString(char *strg)
{
  ch = new char[strlen(strg)+1];

  sprintf(ch,strg);

  if (debug) std::cout << " MyString constructor\n\n";
  
  //std::cout << ++nObj << "\n\n";

  return;
}



MyString::~MyString()
{
  if (ch != NULL) delete [] ch; ch = NULL;

  if (debug) std::cout << " MyString destructor\n\n";
  
  //std::cout << --nObj << "\n\n";

  return;
}



bool MyString::operator==(MyString &strg2)
{
  if (ch == NULL || strg2.ch == NULL) return false;
	
  if (strcmp(ch,strg2.ch)==0) return true;

  return false;
}



bool MyString::operator!=(MyString &strg2)
{
  if (ch == NULL || strg2.ch == NULL) return true;
  
  if (strcmp(ch,strg2.ch)==0) return false;

  return true;
}



bool MyString::operator==(char *strg2)
{
  if (ch == NULL || strg2 == NULL) return false;

  if (strcmp(ch,strg2)==0) return true;

  return false;
}



bool MyString::operator!=(char *strg2)
{
  if (ch == NULL || strg2 == NULL) return true;

  if (strcmp(ch,strg2)==0) return false;

  return true;
}



bool MyString::begins(char *strg2)
{
   if (ch == NULL) return false;

   int l = strlen(strg2);

   if (l > length()) return false;
  
   if (strncmp(ch,strg2,l) == 0) return true;

   return false;   
}



bool MyString::begins(MyString &strg2)
{
   return begins(strg2.asCharArray()); 
}



int MyString::containsAtPos(char *word)
{
   if (ch == NULL || strlen(ch) < 1) return -1;

   if (word == NULL || strlen(word) < 1) return -1;

   int i = 0, l1 = strlen(ch), l2 = strlen(word);

   while (i <= l1-l2 && strncmp(&ch[i],word,l2) != 0) i++;
   
   if (i <= l1-l2) return i;

   return -1;
}



int MyString::containsAtPos(MyString &word)
{
   return containsAtPos(word.ch);
}



bool MyString::contains(char *word)
{
   if (containsAtPos(word) > -1) return true;

   return false;
}



bool MyString::contains(MyString &word)
{
   return contains(word.ch);
}



MyString &MyString::append(MyString &strg2)
{
  return append(strg2.ch);
}



MyString &MyString::append(char *strg2)
{
  if (strg2 == NULL)  return *this;

  int i, l1 = length(), l2 = strlen(strg2);
  if (ch == NULL) l1 = 0;
  
  char *tmp = new char[l1 + l2 + 1];
  
  for (i=0; i<l1; i++) tmp[i]    = ch[i]; 
  for (i=0; i<l2; i++) tmp[i+l1] = strg2[i];  
  tmp[l1+l2] = '\0';

  if (ch != NULL) delete [] ch;

  ch = tmp;

  return *this;  
}



MyString &MyString::append(char c)
{
  char tmp[2]; tmp[0] = c; tmp[1] = '\0';  return append(tmp);
}




MyString &MyString::insert(int pos, char *word)
{
  if (ch == NULL && pos > 0) 
    { std::cout << " WARNING! 'MyString::insert failed!\n\n"; return *this; }
  
  if (word == NULL || strlen(word) < 1) return *this;

  int i, lw = strlen(word), l = strlen(ch);

  if (pos > l) 
    { std::cout << " WARNING! 'MyString::insert failed!\n\n"; return *this; }
	
  char *tmp = new char[lw + l + 1];

  for (i=0; i<pos; i++) tmp[i]     = ch[i];
  for (i=0; i<lw ; i++) tmp[i+pos] = word[i];
  for (i=pos; i<=l;i++) tmp[i+lw]  = ch[i];

  delete [] ch; ch = tmp;
	
  return *this;
}



MyString &MyString::insert(int pos, MyString &word)
{ 
  return insert(pos,word.asCharArray());
}



MyString &MyString::trunc(int i)
{
  if (ch == NULL) return *this;

  if (length() <= i) return *this;

  char *tmp = new char[i+1]; 

  for (int j=0; j<i; j++) tmp[j] = ch[j]; tmp[i] = '\0';

  delete [] ch; ch = tmp;

  return *this;
}



MyString &MyString::free(void)
{
  if (ch != NULL) delete [] ch; 
  
  ch = NULL;

  return *this;
}



char *MyString::asCharArray(void)
{
   return ch;
}


char &MyString::operator[](int i)
{
  if (ch == NULL || strlen(ch)-1<i) 
    { std::cout << " ERROR! MyString operator [],  i > length\n\n"; exit(0); }
  
  return ch[i];
}



int MyString::length(void)
{
   if (ch == NULL) return -1;
   
   return strlen(ch);
}



MyString &MyString::read(std::istream &stream)
{
  if (ch != NULL) { delete [] ch; ch = NULL; }

  if (!stream) return *this;
  
  char c, line[MAX_READ_LINE_LENGTH];
      
  int i = 0;

  while (1)
  {
    while (stream.get(c) && c != APPEND_NEXT_LINE_MARK && c != '\n' && i<MAX_READ_LINE_LENGTH-1)  
      line[i++] = c;

    if (c == APPEND_NEXT_LINE_MARK) while (stream.get(c) && c != '\n'); else break;
  }
  
  line[i] = '\0';
  
  if (c!='\n') 
  {
     if (debug) std::cout << " reading a very long line ... " << line << "  <" << c << ">\n\n";
     
     this->append(line).append(c);

     while (1)
     {
       while (stream.get(c) && c != APPEND_NEXT_LINE_MARK && c != '\n')  this->append(c);

       if (c == APPEND_NEXT_LINE_MARK) while (stream.get(c) && c != '\n'); else break;
     }

     return *this;
  }
 
  int l = strlen(line) + 1;

  ch = new char[l];

  for (i=0;i<l;i++) ch[i] = line[i]; 

  return *this;	
}
   


MyString &MyString::getNextLine(std::istream &stream)
{
  while (stream && read(stream).strip().countWords() < 1);

  return *this;
}



MyString &MyString::input(void)
{
  return read(std::cin);
}



MyString &MyString::inputKeepIfReturn(void)
{
  std::string tmp;

  getline(std::cin,tmp);
   
  int l = tmp.length();

  if (l == 0) return *this;
  
  if (ch != NULL) delete [] ch;
  
  ch = new char[l+1];

  for (int i=0;i<l;i++) ch[i] = tmp[i]; ch[l] = '\0';
 
  return *this;	
}



int MyString::which(char **list)
{
  if (ch == NULL) return -1;

  int i = 0;   while (list[i] != NULL && strcmp(ch,list[i]) != 0) i++;

  if (list[i] == NULL) return -1;
  
  return i;
}



int MyString::which(MyString *list)
{
  int i = 0;
  
  while (list[i].ch != NULL && *this != list[i]) i++;

  if (list[i].ch == NULL) return -1;
  
  return i;
}



void MyString::print(void)
{
  if (ch != NULL)  std::cout << ch;
}



#ifndef WINDOWS
  #define DIRECTORY_LETTERS "/1234567890ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_~."
#else
  #define DIRECTORY_LETTERS "\\1234567890ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_:."
#endif
#define FILE_LETTERS "1234567890ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_."

bool MyString::FileOK(void)
{
  if (ch == NULL) return false;
	
  int i = 0, l = length();
  
  if (l<1) return false;

  while (i<l && charInCharArray(ch[i],FILE_LETTERS)) i++;

  if (i < l) return false;
  
  return true;
}




bool MyString::DirectoryOK(void)
{
  if (ch == NULL) return false;
	
  int i = 0, l = length();

  if (l<1) return false;

  while (i<l && charInCharArray(ch[i],DIRECTORY_LETTERS)) i++;

  if (i < l) return false;
  
  return true;
}



MyString &MyString::strip(void)  // remove comments, remove leading and trailing spaces
{
  if (ch == NULL) return *this;

  int i = 0, l = length();
  
  while (i < l && !charInCharArray(ch[i],COMMENT_MARK)) i++;

  trunc(i);

  i = strlen(ch); while (ch[i-1] == ' ') i--; trunc(i);
  
  i = 0; while (ch[i] == ' ') i++; 
  if (i>0)
    {  char *tmp = new char[strlen(ch)-i+1]; 
       for (int j=0; j<strlen(ch)-i+1; j++) tmp[j] = ch[i+j]; 
       delete [] ch;  ch = tmp; }

  if (strlen(ch) < 1) { delete [] ch; ch = NULL; }
  
  return *this;
}



MyString &MyString::stripToMin(void) // remove multiple WEAK_WORD_SEPARATORS
{
  strip();
 
  if (ch == NULL) return *this;
  
  char *tmp;
  int first = strlen(ch), j, l, i = 0; 
  while (i<strlen(ch))
  {
     if (charInCharArray(ch[i],ALL_WORD_SEPARATORS)) first = ++i;
     while (charInCharArray(ch[i],WEAK_WORD_SEPARATORS)) i++;
     if (charInCharArray(ch[i],ALL_WORD_SEPARATORS) && 
	 charInCharArray(ch[first-1],WEAK_WORD_SEPARATORS)) first--;
     if (i > first)
     {
	l = strlen(ch) - i + first;
	//std::cout << first << " " << i << " " << l << "\n\n";
        tmp = new char[l+1];
	for (j=0; j<first; j++) tmp[j] = ch[j];
	for (j=first; j<l+1; j++) tmp[j] = ch[j+i-first];
	delete [] ch; 
	ch = tmp; 
        i = first-1;
     }
     first = strlen(ch)+5;
     i++;
  }

  return *this;  
}



int MyString::countWords()
{
  if (ch == NULL || strlen(ch) < 1) return 0;
  
  int i = 0, n = 0, l = strlen(ch), m;

  while (i <= l) 
    {  
       while (charInCharArray(ch[i],WEAK_WORD_SEPARATORS)) i++;
       //std::cout << ch[i] << i << " first\n";

       while (!charInCharArray(ch[i],ALL_WORD_SEPARATORS) && ch[i] != '\0') 
       {
         if (ch[i] == '{')
         {
           m = i++;
           while (ch[i] != '}' && ch[i] != '\0') i++;
           if (ch[i] == '\0') i = m;
         }
         i++;
       }
      
       //std::cout << ch[i] << i << " last\n";

       n++;

       while (charInCharArray(ch[i],WEAK_WORD_SEPARATORS)) i++;

       if (charInCharArray(ch[i],ALL_WORD_SEPARATORS) || ch[i] == '\0') i++; 
    }
 
  return n;	
}




bool MyString::copyAfter(char c, MyString &word2)
{ 
  if (ch == NULL) return false;
	
  int i = 0;
  
  while (ch[i] != c && ch[i] != '\0') i++;

  if (ch[i] == '\0') return false;

  if (ch[i+1] == '\0') return false;
  
  word2.free().append(&(ch[i+1]));

  return true;
}




int MyString::split(char ***words)
{
  if (ch == NULL) return 0;

  if (length() == 0) return 0;
 
  int i = 0, n = countWords(), first, l = strlen(ch), j, lw, m;
  
  *words = new char* [n+1];  // allocate memory for pointer array

  i = 0, n = 0;
  
  while (i <= l) 
    {  
       while (charInCharArray(ch[i],WEAK_WORD_SEPARATORS)) i++;
       //std::cout << ch[i] << i << " first\n";
       first = i;

       while (!charInCharArray(ch[i],ALL_WORD_SEPARATORS) && ch[i] != '\0')
       {
         if (ch[i] == '{')
         {
           m = i++;
           while (ch[i] != '}' && ch[i] != '\0') i++;
           if (ch[i] == '\0') i = m;
         }
         i++;
       }
       //std::cout << ch[i] << i << " last\n";
       
       lw = i - first;

       (*words)[n] = new char [lw + 1];

       for (j=0; j<lw; j++) (*words)[n][j] = ch[first + j];  (*words)[n][lw] = '\0';
       
       n++;

       while (charInCharArray(ch[i],WEAK_WORD_SEPARATORS)) i++;

       if (charInCharArray(ch[i],ALL_WORD_SEPARATORS) || ch[i] == '\0') i++; 
    }
  
  (*words)[n] = NULL;

  if (debug) std::cout << " Don't forget to free the memory allocated here for 'words'!\n\n";
      	
  return n;
}



int MyString::split(MyString **words)
{
   char **tmp;

   int i, n = split(&tmp);

   if (n == 0) return 0;
   
   *words = MyStringArray(tmp);

   return n;
}



bool MyString::toInt(int *i, bool errMsg) 
{
  if (!scanInt(ch,i,errMsg)) return false;

  return true;
}



bool MyString::toDbl(double *x, bool errMsg) 
{
  if (!scanDbl(ch,x,errMsg)) return false;

  return true;
}



bool MyString::toInt(myVector<int> &a, bool errMsg) 
{
  a.free();

  MyString *word;

  char *ch0 = ch, *chi;

  int l = strlen(ch);

  if (ch[0] == '{') ch++;
  if (ch0[l-1] == '}') ch0[l-1] = '\0';

  int nw = this->split(&word), i = 0, j, x;

  while (i < nw)
  {
    j = word[i].containsAtPos("-");

    chi = word[i].ch;

    if (j > 0)
    { 
      chi[j] = '\0';

      if (!scanInt(chi,a.append(),errMsg)) goto error;

      if (!scanInt(chi+j+1,&x,errMsg)) goto error;

      for (j=a.lastCoeff()+1; j<x; j++) a.append(j);

      a.append(x);
    }
    else if (!scanInt(chi,a.append(),errMsg)) goto error;

    i++;
  }
  delete [] word; ch = ch0;

  return true;

  error:

  delete [] word; ch = ch0;

  return false;
}





// overload printing to screen

std::ostream &operator<<(std::ostream &os,MyString &strg)
{
   if (strg.asCharArray() != NULL) os << strg.asCharArray(); return os;
}



// overload writing to file

std::ofstream &operator<<(std::ofstream &os,MyString &strg)
{
   if (strg.asCharArray() != NULL) os << strg.asCharArray(); return os;
}



// generate an array of MyStrings from  char **list

MyString *MyStringArray(char **list)
{
  int i = 0;
	
  int c = 0; while (list[c] != NULL) c++;

  MyString *tmp = new MyString[c];
 
  for (i=0; i<c; i++)  { tmp[i].append(list[i]); delete [] list[i]; } delete [] list;
  
  return tmp;
}




