
#ifndef incl_FunctionsProgram_h
#define incl_FunctionsProgram_h

#include <fstream>
#include <iostream>

#include "MyStringList.h"
#include "Macro.h"
#include "MacroCommand.h"
#include "MathVector.h"



// functions in program

void       play                          ();
void       prgPrintLogo                  ();
bool       prgFileExist                  (char*, char*);
bool       prgFileExist                  (char*);
bool       prgFileExist                  (MyString &, MyString &);
bool       prgFileExist                  (MyString &);
bool       prgReadFile                   ();
bool       prgReadFileNames              ();
//bool       prgLastProject                ();
//void       prgWriteLastProject           ();

inline  void  prgError(int id, char* fname, char* msg)
{
  std::cout << " ERROR (" << id << ") in " << fname << ": " << msg << "\n\n";

  exit(0);
}


inline  void  prgWarning(int id, char* fname, char* msg)
{
  std::cout << " WARNING (" << id << ") in " << fname << ": " << msg << "\n\n";

  return;
}

void       prgLoadMacros                 ();
int        prgNextDomainKey              (std::ifstream &);
bool       prgNextDataKey                (std::ifstream &, MyString &);
int        prgListCharArrayLength        (char**);
void       prgResolveIncludeFile         (std::ofstream &, char*);
void       prgReadRunControl             (std::ifstream &);
void       prgReadTimeFunctions          (std::ifstream &);
bool       prgStringToMacroCmd           (MacroCommand &, MyString&);
void       prgMacroCmdToSimpleString     (MyString&, MacroCommand &);
void       prgGetMacroCommand            (MacroCommand &, int &);
void       prgExecNextRunCtrlCmd         ();
void       prgCloseProject               ();

//void       prgPrintSimpleMatrix          (double*, int, int, int, int, bool, int indent = 0, bool numbers = false);
//void       prgPrintSimpleMatrixDiff      (double*, double*, int, int, int, int, bool, int indent = 0, bool numbers = false);
//void       prgPrintSimpleMatrixSec       (double*, int, int, int, int, int, int, int, bool, int indent = 0, bool numbers = false);
//void       prgPrintSimpleMatrixSecDiff   (double*, double*, int, int, int, int, int, int, int, bool, int indent = 0, bool numbers = false);
//void       prgCompareTwoSimpleMatrices   (double*, double*, char*, char*, char*, int, int, int, int, bool, int indent = 0, bool interactive = false, bool numbers = false);
//char       prgUpdateDisplay              (bool wait = false);

inline  bool  prgNAN(double x)
{
  char tmp[20];
  double y;

  sprintf(tmp,"%12.5g  ",x);

  while(tmp[0] == ' ')
  {
    for(int i=0; i<strlen(tmp); i++)
	  tmp[i] = tmp[i+1];
  }

  while(tmp[strlen(tmp)-1] == ' ')
  {	tmp[strlen(tmp)-1] = '\0';  }

  if(!scanDbl(tmp,&y,false))
    return true;

  return false;
}



inline  bool  prgNAN(double *x, int n)
{
  for(int i=0; i<n; i++) if (prgNAN(x[i])) return true;

  return false;
}


void       prgMacroTest                  (double, double);
bool       prgReadColonSepListVectorInt  (std::ifstream &, MyString &, List< Vector<int> > &);
bool       prgReadLnBrkSepListVectorInt  (std::ifstream &, MyString &, List< Vector<int> > &);
bool       prgReadLnBrkSepListVectorDbl  (std::ifstream &, MyString &, List< Vector<double> > &);
//void       prgTest                       (void);
void       prgCheckCompressedMatrix      (int, int, VectorArray<double> &, VectorArray<int> &, VectorArray<int> &, VectorArray<int> &);
void       prgPlotMatrixPattern          (int, int, VectorBase<int> &, VectorBase<int> &, char* fileName = NULL);
int        prgWriteToTFile               (MyString &);
int        prgWriteToTFile               (char *);




// for Fortran subroutines
extern "C" 
{
  void prgerror_(int *id, char* fname, char* msg);

  void prgwarning_(int *id, char* fname, char* msg);
}


#endif

