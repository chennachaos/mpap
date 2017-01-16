
#ifndef incl_MacroCommand_h
#define incl_MacroCommand_h


#include "MyString.h"
#include "MathVector.h"
#include "Parameter.h"
#include "List.h"


class MacroCommand: public ListItem
{
  public:
	  
    int            ii;

    MyString       strg;

    ListInfinite<Parameter> p;
    
    bool operator!(void) { if (ii < 0) return true; return false; }

    void free(void) { ii = -1; strg.free(); p.free(); return; }
};



#endif

