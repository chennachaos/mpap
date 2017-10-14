
#ifndef incl_MyStringList_h
#define incl_MyStringList_h

#include <iostream> 

#include "List.h" 
#include "MyString.h" 


class MyStringList: public List<MyString>
{
  public:
    
    MyStringList(void);
    virtual ~MyStringList();

    void addNew(char*, char* s1 =NULL, char* s2 =NULL, char* s3 =NULL, char* s4 =NULL, 
	               char* s5 =NULL, char* s6 =NULL, char* s7 =NULL, char* s8 =NULL, 
		       char* s9 =NULL, char* s10=NULL, char* s11=NULL, char* s12=NULL, 
		       char* s13=NULL, char* s14=NULL, char* s15=NULL, char* s16=NULL, 
		       char* s17=NULL, char* s18=NULL, char* s19=NULL, char* s20=NULL, 
		       char* s21=NULL, char* s22=NULL, char* s23=NULL, char* s24=NULL, 
		       char* s25=NULL, char* s26=NULL, char* s27=NULL, char* s28=NULL, 
		       char* s29=NULL, char* s30=NULL, char* s31=NULL, char* s32=NULL, 
		       char* s33=NULL, char* s34=NULL, char* s35=NULL, char* s36=NULL, 
		       char* s37=NULL, char* s38=NULL, char* s39=NULL, char* s40=NULL ); 
   
    void print();
   
    char **genCharArrays(void);

    int which(char *);
    int which(MyString &);
    
    int whichBegins(char *);
    int whichBegins(MyString &);
    
   private:

};


std::ostream &operator<<(std::ostream  &, MyStringList &);



#endif
