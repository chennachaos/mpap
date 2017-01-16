
#include <Xm/Xm.h>

#include "FunctionsProgram.h"
#include "FunctionsUnix.h"


void grpOpen(Widget w, XtPointer client_data, XtPointer call_data)
{
  Widget dummy = NULL;

/*
  for(int ii=0;ii<20;ii++)
  {
    cout << " AAAAAAAAAAAAAAAAA " << endl;
    cout << endl;
  }
*/
  
  if (prgLastProject()) grpPopupLastProject(dummy, NULL, NULL);
  
    else  grpPopupFileSelect (dummy, NULL, NULL);

  return;
}


