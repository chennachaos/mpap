
#ifndef incl_Macro_h
#define incl_Macro_h


#include "List.h"
#include "MyString.h"
#include "MacroDialogData.h"
#include "MathVector.h"


class Macro: public ListItem
{ 
  public:

    Macro(void);
 
    Macro(int);                      // constructor

    virtual ~Macro();                // destructor
    
    MyString     type;               // macro type
    MyString     name;               // macro name
    MyString     what;               // short description (one line)
    VectorFixed<bool,10> sensitivity;// sensitivity flags

    bool operator!(void);            // if macro not defined ...
    
    void       init(void);           // initialise macro
    int        exec(void);           // execute macro

    MacroDialogData db;              // dialog box data

    int        id;                   // macro id to identify function macro1,2,3,..

    VectorInfinite<double> p;
    MyString               strg;

  private:                                

};




inline bool Macro::operator!(void)
{
   if (!name) return true;

   return false;
}



std::ostream  &operator<<(std::ostream  &, Macro &);




#include <iostream>         // Doing this here,
#include "FunctionsMacro.h" // I avoids typing this 
#include "RunModeEnum.h"    // in all the macro functions
#include "MathBasic.h"

#endif




