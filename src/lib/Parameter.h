
#ifndef incl_Parameter_h
#define incl_Parameter_h


#include "MyString.h"
#include "List.h"


enum NumberType { NON_NEG_INT, NON_NEG_REAL, ALL_INT, ALL_REAL };



class Parameter: public ListItem
{
  public:

    Parameter(void);

    ~Parameter(void);
	  
    Parameter &operator=(Parameter &);
    
    void free(void);
 
    bool interprete(char*, NumberType);
    
    double evaluate(int);
    
    double    x;
	  
    MyString  fctName, fctPar;
};





#endif

