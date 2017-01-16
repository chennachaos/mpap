
#ifndef incl_UnixFunctions_h
#define incl_UnixFunctions_h

#include <Xm/Xm.h>




// callback functions for main user interface

void   grpOpen                     (Widget, XtPointer, XtPointer);
void   grpClose                    (Widget, XtPointer, XtPointer);
void   grpExit                     (Widget, XtPointer, XtPointer);
void   grpPopupHelp                (Widget, XtPointer, XtPointer);
void   grpPopupAbout               (Widget, XtPointer, XtPointer);
void   grpPopupMacro               (Widget, XtPointer, XtPointer);
void   grpDrawingAreaExposed       (Widget, XtPointer, XtPointer);
void   grpDrawingAreaResized       (Widget, XtPointer, XtPointer);
void   grpCommandListToText        (Widget, XtPointer, XtPointer);
void   grpCommandEntered           (Widget, XtPointer, XtPointer);
void   grpDrawingAreaMouseInput    (Widget, XtPointer, XtPointer);

                                  
// other callbacks                
                                  
void   grpPopdownMacro             (Widget, XtPointer, XtPointer);
void   grpMacroHelp                (Widget, XtPointer, XtPointer);
                                  
void   grpPopupLastProject         (Widget, XtPointer, XtPointer);
void   grpPopdownLastProject       (Widget, XtPointer, XtPointer);
                                  
void   grpPopupFileSelect          (Widget, XtPointer, XtPointer);
void   grpPopdownFileSelect        (Widget, XtPointer, XtPointer);


// other

void   grpOpenProject              (void);

// general help functions

Widget grpGenThreeButtonControlArea(Widget*, Widget*, Widget*, Widget, char*, char*, char*);
Widget grpGenTwoButtonControlArea  (Widget*, Widget*, Widget, char*, char*);
Widget grpFindManagerWidget        (Widget);
void   grpPopdownDialogShell       (Widget);




#endif

