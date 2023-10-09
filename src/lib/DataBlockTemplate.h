
#ifndef incl_DataBlockTemplate_h
#define incl_DataBlockTemplate_h

#include "MyStringList.h"
#include "MathVector.h"


class DataBlockTemplate
{
  public:

    DataBlockTemplate(void);

    ~DataBlockTemplate();
   
    bool initialise(char *,     char **l1=NULL, char **l2=NULL, char **l3=NULL, char **l4=NULL);
    bool initialise(MyString &, char **l1=NULL, char **l2=NULL, char **l3=NULL, char **l4=NULL);
    
    bool expandToMatch(DataBlockTemplate &);

    bool readBlock(std::ifstream &, MyString &, 
                   myVector<int> &, myVector<double> &, MyStringList &, myVector<int> &);

    void free(void);
    
  private:
    
    bool readLine (std::ifstream &, MyString &, 
                   myVector<int> &, myVector<double> &, MyStringList &, myVector<int> &);
    
    myVector<int>   type, expr, t, c;

    myVector<int>    iDef, lDef;
    myVector<double> dDef;
    MyStringList   sDef;

    char*** list;

    int  cnt;
    
    void append(int, int);
    
};


#endif

