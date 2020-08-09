
#ifndef incl_MacroList_h
#define incl_MacroList_h


#include "Macro.h"
#include "List.h"



class MacroList: public List<Macro>
{ 
  public:

    MacroList(void);       // constructor
    virtual ~MacroList();  // destructor
    
    int ntype;             // number of macro types
    int *jp;               // profile of list for macro types

    int checkNames(void);  // make sure macro names and types are distinct
    void sort(void);       // sort macro list
                                          
  private:                                

};



#endif




