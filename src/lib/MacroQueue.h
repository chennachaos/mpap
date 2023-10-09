
#ifndef incl_MacroQueue_h
#define incl_MacroQueue_h


#include "List.h"
#include "MacroCommand.h"
#include "MathVector.h"
#include "Loop.h"
#include "If.h"



class MacroQueue
{
  public:

    MacroQueue(void);
	  
    void append(MacroCommand &);

    void reset(void);

    List<Loop> loop;

    List<If> iff;
    
    List<MacroCommand> macCmd;

    myVector<double> p;      // double parameters of the last macro to be called
    MyString       strg;   // string parameter  of the last macro to be called
                           // (to be set in the macro and picked up in grpDrawingAreaMouseInput)
                              
  private:
   
    int openLoops, openIfs;

};

#endif


