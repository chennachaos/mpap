
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
                   Vector<int> &, Vector<double> &, MyStringList &, Vector<int> &);

    void free(void);
    
  private:
    
    bool readLine (std::ifstream &, MyString &, 
                   Vector<int> &, Vector<double> &, MyStringList &, Vector<int> &);
    
    Vector<int>   type, expr, t, c;

    Vector<int>    iDef, lDef;
    Vector<double> dDef;
    MyStringList   sDef;

    char*** list;

    int  cnt;
    
    void append(int, int);
    
};


#endif

