
#ifndef incl_EssGrpFunctions_h
#define incl_EssGrpFunctions_h


#include "MyString.h"


// functions which must be provided with exactly this interface 
// under both linux and windows


void essGrpStartGUI                     (int, char**);
void essGrpWriteTime                    (void);
void essGrpWriteStatus                  (void);
void essGrpWriteProject                 (char *strg = NULL);
void essGrpSetMacroSens                 (void);
void essGrpUpdateDisplay                (void);
void essGrpGetPlotAreaGeom              (void);
void essGrpDrawLine                     (int, int, int, int);
void essGrpWipe                         (void);
void essGrpSetColour                    (int, bool);
void essGrpSetSensAllButDrawingArea     (bool);
void essGrpFillPoly                     (int*, int);
void essGrpFillCircle                   (int*, int*);
void essGrpDrawCircle                   (int*, int*);
void essGrpPutText                      (int*, char*, int, bool BGFlag = false);
void essGrpCheckForResizeEvents         (void);
void essGrpSetCursor                    (bool);
int  essGrpWaitForMouseButtonPressed    (void);
void essGrpCopyPixmap                   (void);
bool essGrpInquire                      (char*, char*, char*, MyString*);
int  essGrpWasMouseButtonPressed        (int *x = NULL, int *y = NULL, bool *inDrawingArea = NULL);


#endif

