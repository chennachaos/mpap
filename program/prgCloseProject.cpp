//
#include "DomainTree.h"
#include "Files.h"
#include "RunControl.h"
#include "MacroQueue.h"
#include "MpapTime.h"
////#include "Plot.h"
#include "List.h"
#include "TimeFunction.h"
#include "ComputerTime.h"
#include "Counter.h"
#include "PlotVTK.h"


extern DomainTree            domain;
extern Files                 files;
extern RunControl            runCtrl;
extern MacroQueue            macroQueue;
extern MpapTime              mpapTime;
//extern Plot                  plot;
extern List<TimeFunction>    timeFunction;
extern ComputerTime          computerTime;
extern ListInfinite<Counter> counter;
extern PlotVTK  plotvtk;

using namespace std;


void prgCloseProject(void)
{
  if (runCtrl.mode == NOPROJ) return;
      
  mpapTime.reset();

  files.reset();
  
  macroQueue.reset();

  //plot.reset();
  
  runCtrl.reset();   
  
  runCtrl.newMode(NOPROJ);
  runCtrl.newStatus(NOPROJECT);
  
  //essGrpWriteTime();
  //essGrpWriteProject();

  timeFunction.free();

  domain.free();

  domain.reset();

  computerTime.reset();

  counter.free();
  
  plotvtk.reset();

  return;
}


