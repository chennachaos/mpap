
#ifndef incl_RunControl_h
#define incl_RunControl_h

#include "MyStringList.h"
#include "List.h"
#include "RunModeEnum.h"
#include "MyStatusEnum.h"
//#include "FunctionsEssGrp.h"
#include "MathVector.h"
////#include "Plot.h"


////extern Plot plot;


struct RunControl

{ 
  RunControl(void) { reset(); return; }
  
  int act;                      // active command currently being executed

  Vector<int> cmd;              // array of integer command keys

  List<MyStringList> batch;     // List of of Batch macro lists

  RunMode mode;                 // current run mode

  MyStatus status, fixedStatus; // current status

  int macCnt;                   // macro counter

  bool quit;                    // exit flag for use with GUI
 

  // functions to control the status

  void reset(void) { act = -1; 
	             macCnt = 0; 
		     batch.free(); 
		     cmd.free(); 
		     quit = false; 
	             fixedStatus = UNDEF; 
		     return; }
  
  void newMode(RunMode rm)
  {
    mode = rm;
    if (mode == BATCH) 
       cout <<  " BATCH mode " << endl;
        //plot.suppressCopyPixmap = true;
    else 
    { 
        //plot.suppressCopyPixmap = false;
        //essGrpCopyPixmap();
    }
    //essGrpSetMacroSens();
    //if (rm == PRE)
    //essGrpWriteProject("pre-processor");

    return;
  }

  void newStatus(MyStatus st)
  {
    status = st;
    //essGrpWriteStatus();
    return; 
  }

  void fixStatus(MyStatus st)
  {
    fixedStatus = st;
    //essGrpWriteStatus();
    return;
  }
  
  void freeFixedStatus()
  {
    fixedStatus = UNDEF;
    //essGrpWriteStatus();
    return;
  }
  
  void busy(void)
  {
    //essGrpSetCursor(true);
    if (fixedStatus != UNDEF)
      return;
    if (mode == INTER)
      newStatus(INTERACTIVEBUSY);
    return;
  }

  void notBusy(void)
  {
    //essGrpSetCursor(false);
    if (fixedStatus != UNDEF)
      return;
    if (mode == PRE || mode == INTER)
      newStatus(INTERACTIVE); 
    if (mode == NOPROJ)
      newStatus(NOPROJECT); 

    return; 
  }

};

#endif

