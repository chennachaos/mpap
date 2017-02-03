
#ifndef incl_Global_h
#define incl_Global_h


#include "Files.h"
#include "DomainTree.h"
#include "MacroList.h"
#include "RunControl.h"
#include "MpapTime.h"
#include "MacroQueue.h"
#include "WorkSpace.h"
#include "SolverWorkSpace.h"
#include "List.h"
#include "TimeFunction.h"
#include "ComputerTime.h"
#include "Counter.h"
#include "SolverTime.h"

DomainTree domain;

MacroList  macro;

bool noGUI           = false;
bool batch           = false;
bool keep            = false;
bool aspectRatioCorr = false;
bool lastProj        = false;
bool test            = false;

Files files;

RunControl runCtrl;

MacroQueue macroQueue;

MpapTime mpapTime;

ComputerTime computerTime;

SolverTime      solverTime;
SolverWorkSpace solverWorkSpace;
WorkSpace       workSpace;

List<TimeFunction> timeFunction;

ListInfinite<Counter> counter;

double globalMaxIncrement;

void *macro2mousePtr;

int countThis = 0;

//int *mpapVar;

// for debugging

bool debug = false;
bool wulf  = false;
bool sony  = false;
bool deniz = false;
bool readIfileInfo = false;
bool vtkFlag=false;

// stuff for debugging

int  nObj = 0;
int  nVector = 0;
bool printThisRubbish = false;

#endif



