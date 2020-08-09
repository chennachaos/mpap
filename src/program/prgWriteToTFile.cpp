
#include "MyString.h"
#include "MpapTime.h"
#include "TimeFunction.h"
#include "Files.h"
#include "Definitions.h"
#include "List.h"
//#include "Global.h"


extern MpapTime           mpapTime;
extern List<TimeFunction> timeFunction;
extern Files              files;



int prgWriteToTFile(MyString &addStr)
{
  char fct[] = "prgWriteToTFile";

  int         i;
  char        tmp[20];
  MyString    timeStr, fname;

  if (mpapTime.write + 1.e-12 < mpapTime.cur) 
  { 
    sprintf(tmp,"\n%12.5f",mpapTime.cur);
    timeStr.append(tmp);
    for (i=0; i<timeFunction.n; i++)
    {
      sprintf(tmp," %12.5f",timeFunction[i].prop);
      timeStr.append(tmp);
    }
    mpapTime.write = mpapTime.cur;

    // open Tfile

    if (!files.Tout.is_open()) 
    {
      fname.append(files.projDir).append(SLASH).append(files.Tfile);
      files.Tout.open(fname.asCharArray(),ios_base::app);
    }

    // write to Tfile

    files.Tout << ' ' << timeStr.stripToMin() << ' ' << addStr.stripToMin();
  }
  else
  {
    // open Tfile

    if (!files.Tout.is_open()) 
    {
      fname.append(files.projDir).append(SLASH).append(files.Tfile);
      files.Tout.open(fname.asCharArray(),ios_base::app);
    }

    // write to Tfile

    files.Tout << ' ' << addStr.stripToMin();
  }

  // close Tfile

  files.Tout.close();
  
  return 0;
}







int prgWriteToTFile(char *addStr)
{
  MyString tmp;

  tmp.append(addStr);

  return prgWriteToTFile(tmp);
}



