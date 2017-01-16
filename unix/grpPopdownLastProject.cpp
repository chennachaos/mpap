
#include <Xm/Xm.h>

#include "FunctionsProgram.h"
#include "FunctionsUnix.h"
#include "Files.h"


extern Files files;


void grpPopdownLastProject(Widget w, XtPointer client_data, XtPointer call_data)
{
  Widget dummy = NULL;
  
  if (strcmp(XtName(w),"Cancel") == 0) { grpPopdownDialogShell(grpFindManagerWidget(w)); return; }
  
  if (strcmp(XtName(w),"Browse") == 0 || !prgFileExist(files.projDir,files.Ifile)) 
      {
          grpPopdownDialogShell(grpFindManagerWidget(w));
	  grpPopupFileSelect(dummy, NULL, NULL);
	  return;  
      }
  
  grpPopdownDialogShell(grpFindManagerWidget(w));

  grpOpenProject();
 
  return;  
}

 

