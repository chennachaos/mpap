

#include "MacroQueue.h"
#include "MacroList.h"
#include "MathBasic.h"
#include "RunControl.h"
//#include "FunctionsEssGrp.h"
#include "FunctionsProgram.h"
////#include "Plot.h"
#include "Files.h"


extern bool       noGUI;
extern MacroList  macro;
extern RunControl runCtrl;
////extern Plot       plot;
extern Files      files;


using namespace std;


MacroQueue::MacroQueue()
{
   reset();

   return;
}




void MacroQueue::append(MacroCommand &mC)
{
  int i, ii, j, act, rtrn;

  char tmp[30];

  // open new loop
  
  if (macro[mC.ii].name == "loop")
  {
    //cout << "loop\n";

    for (i=mC.p.n; i>0; i--) mC.p[i] = mC.p[i-1];
    
    loop.add(new Loop);
    
    mC.p[0].free();
    i            = loop.n - 1;
    mC.p[0].x    = (double) i;
    
    loop[i].beg  = macCmd.n;
    loop[i].end  = -1;
    loop[i].cnt  = 1;
    loop[i].name = mC.strg;
    
    if (loop[i].name == "")  
    { 
      if (runCtrl.mode == BATCH) j = 0; else j = 1;
      sprintf(tmp,"loop%d",i+j); 
      loop[i].name = tmp; 
    }
   
    j = 0; while (j<i) if (loop[j].name == loop[i].name) 
      { loop.trunc(i); COUT << "ignoring 'loop'! duplicate loop name!\n\n"; return; }
    else j++;

    openLoops++;
  }

  // close loop
  
  else if (macro[mC.ii].name == "next")
  {
    //cout << "next\n";

    i = loop.n - 1;  while (loop[i].end >= 0) i--;
    
    if (openIfs > 0) 
    {
      j = iff.n - 1; while (iff[j].end >=0) j--;
      if (iff[j].beg > loop[i].beg) 
        { COUT << "ignoring 'next'! close if-statement first!\n\n"; return; }
    }
    
    loop[i].end = macCmd.n;

    mC.p[0].free();
    mC.p[0].x   = (double) i;
    
    openLoops--;
  }

  // configure 'xlop'

  else if (macro[mC.ii].name == "xlop")
  {
    //cout << "xlop\n";

    if (openLoops < 1) { COUT << "ignoring 'xlop'! no open loops!\n\n"; return; }
	  
    i = loop.n - 1;  while (loop[i].end >= 0) i--;
   
    mC.p[0].free(); 
    mC.p[0].x = (double) i;
  }

  // open new if-statement
  
  else if (macro[mC.ii].name == "if")
  {
    //cout << "if\n";

    for (i=mC.p.n; i>0; i--) mC.p[i] = mC.p[i-1];

    iff.add(new If);
    
    mC.p[0].free();
    i          = iff.n - 1;
    mC.p[0].x  = (double) i;
    
    iff[i].beg = macCmd.n;
    iff[i].els = -1;
    iff[i].end = -1;
    
    openIfs++;
  }
  
  // 'else'
  
  else if (macro[mC.ii].name == "else")
  {
    //cout << "else\n";

    i = iff.n - 1;  while (iff[i].end >= 0) i--;
   
    if (openLoops > 0)
    {
      j = loop.n - 1; while (loop[j].end >=0) j--;
      if (loop[j].beg > iff[i].beg) 
        { COUT << "ignoring 'else'! close loop first!\n\n"; return; }
    }
      
    if (iff[i].els > 0) { COUT << "ignoring multiple 'else'!\n\n"; return; }
    
    iff[i].els = macCmd.n;

    mC.p[0].free();
    mC.p[0].x  = (double) i;
  }
  
  // close if-statement

  else if (macro[mC.ii].name == "ndif")
  {
    //cout << "ndif\n";

    i = iff.n - 1;  while (iff[i].end >= 0) i--;
   
    if (openLoops > 0) 
    {
      j = loop.n - 1; while (loop[j].end >=0) j--;
      if (loop[j].beg > iff[i].beg) 
        { COUT << "ignoring 'ndif'! close loop first!\n\n"; return; }
    }
    
    iff[i].end = macCmd.n;
    if (iff[i].els < 0) iff[i].els = iff[i].end;

    mC.p[0].free();
    mC.p[0].x = (double) i;
    
    openIfs--;
  }
  
  // end 

  else if (macro[mC.ii].name == "end")
  {
    //cout << "end\n";

    if (openLoops>0 || openIfs>0) cout << "\n WARNING! open loop or if statements ignored!\n\n";
	    	  
    if (!!mC.p[0].fctName)        cout << "\n WARNING! parameter function ignored!\n\n";
    
    reset();
   
    if (roundToInt(mC.p[0].x) == 1) runCtrl.quit = true;
    
    else if (!noGUI && runCtrl.act > runCtrl.cmd.n - 2) 
      {    
        std::cout << "\n          Input file run control finished.";
        std::cout << "\n          To close the project, select 'project' -> 'close'";
        std::cout << "\n          To quit Mpap2, select 'project' -> 'exit'\n\n";

        if (runCtrl.mode == PRE) runCtrl.newStatus(INTERACTIVE);

        else { runCtrl.newMode(INTER); runCtrl.newStatus(INTERACTIVE); }
      }
      else  prgExecNextRunCtrlCmd();

    return;
  }
  
  // append macro to queue
 
  macCmd.add(new MacroCommand);
  
  i = macCmd.n - 1;
 
  macCmd[i].ii   = mC.ii;
  macCmd[i].strg = mC.strg;
  for (int j=0; j<mC.p.n; j++) macCmd[i].p[j] = mC.p[j];
  
  if (openLoops > 0 || openIfs > 0) return;

  // execute
      
  runCtrl.busy();	
   
  act = 0;
  while (act < macCmd.n)
  {
    //cout << act << "\n\n";
	  
    ii = macCmd[act].ii;
    
    macro[ii].strg = macCmd[act].strg;
   
    for(i=0; i<macCmd[act].p.n; i++) macro[ii].p[i] = macCmd[act].p[i].evaluate(act);
   
    //essGrpUpdateDisplay();

    //if (macro[ii].type == "plot" && plot.psOpen) 
      //files.Pout << "\n%% **** " << macro[ii].name << " ****\n";

    rtrn = macro[ii].exec();  
    
    if (rtrn < 0) 
    { 
      cout << " last macro: " << macro[ii].name << "\n\n";
      prgError(1,"MacroQueue::append","negative macro return value!");
    }

    if (rtrn > 0) act = rtrn - 1; else act++;
  
    //essGrpCheckForResizeEvents();

    //essGrpUpdateDisplay();
    
    //if (!plot.suppressCopyPixmap) essGrpCopyPixmap();

    macro[ii].strg.free();
    macro[ii].p.free();
  }

  runCtrl.notBusy();		
 
  reset();
  
  return;
}






void MacroQueue::reset(void)
{
  loop.free();
  iff.free();
  macCmd.free();
  openLoops = 0;
  openIfs   = 0;

  return;
}









