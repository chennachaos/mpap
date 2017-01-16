
#ifndef incl_FunctionsProgram_h
#define incl_FunctionsProgram_h

#include <fstream>
#include <iostream>

#include "MyStringList.h"
#include "Macro.h"
#include "MacroCommand.h"
#include "MathVector.h"



// functions in program

void       play                          (void);

void       prgPrintLogo                  (void);
bool       prgFileExist                  (char*, char*);
bool       prgFileExist                  (char*);
bool       prgFileExist                  (MyString &, MyString &);
bool       prgFileExist                  (MyString &);
bool       prgReadFile                   (void);
bool       prgReadFileNames              (void);
bool       prgLastProject                (void);
void       prgWriteLastProject           (void);
void       prgError                      (int, char*, char*);
void       prgWarning                    (int, char*, char*);
void       prgLoadMacros                 (void);
int        prgNextDomainKey              (std::ifstream &);
bool       prgNextDataKey                (std::ifstream &, MyString &);
int        prgListCharArrayLength        (char**);
void       prgResolveIncludeFile         (std::ofstream &, char*);
void       prgReadRunControl             (std::ifstream &);
void       prgReadTimeFunctions          (std::ifstream &);
bool       prgStringToMacroCmd           (MacroCommand &, MyString&);
void       prgMacroCmdToSimpleString     (MyString&, MacroCommand &);
void       prgGetMacroCommand            (MacroCommand &, int &);
void       prgExecNextRunCtrlCmd         (void);
void       prgCloseProject               (void);

void       prgPrintSimpleMatrix          (double*, int, int, int, int, bool, 
		                          int indent = 0, 
				          bool numbers = false);
void       prgPrintSimpleMatrixDiff      (double*, double*, int, int, int, int, bool, 
		                          int indent = 0, 
					  bool numbers = false);
void       prgPrintSimpleMatrixSec       (double*, int, int, int, int, int, int, int, bool, 
		                          int indent = 0, 
				          bool numbers = false);
void       prgPrintSimpleMatrixSecDiff   (double*, double*, int, int, int, int, int, int, int, 
                                          bool, int indent = 0, 
					  bool numbers = false);
void       prgCompareTwoSimpleMatrices   (double*, double*, char*, char*, char*, int, 
		                          int, int, int, bool, 
		                          int indent = 0, 
				          bool interactive = false, 
				          bool numbers = false);
char       prgUpdateDisplay              (bool wait = false);
bool       prgNAN                        (double);
bool       prgNAN                        (double*, int);
void       prgMacroTest                  (double, double);
bool       prgReadColonSepListVectorInt  (std::ifstream &, MyString &, List< Vector<int> > &);
bool       prgReadLnBrkSepListVectorInt  (std::ifstream &, MyString &, List< Vector<int> > &);
bool       prgReadLnBrkSepListVectorDbl  (std::ifstream &, MyString &, List< Vector<double> > &);
void       prgTest                       (void);
void       prgCheckCompressedMatrix      (int, int,
                                          VectorArray<double> &,
                                          VectorArray<int> &,
                                          VectorArray<int> &,
                                          VectorArray<int> &);
void       prgPlotMatrixPattern          (int, int,
                                          VectorBase<int> &, VectorBase<int> &,
                                          char* fileName = NULL);
int        prgWriteToTFile               (MyString &);
int        prgWriteToTFile               (char *);


extern "C" 
{
  void prgerror_       (int*, char*, char*); 
  void prgwarning_     (int*, char*, char*); 
}


#endif

