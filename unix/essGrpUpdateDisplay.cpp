
#include <Xm/Xm.h>

#include "UnixGUI.h"


extern UnixGUI unixGUI;
extern bool    noGUI;


void essGrpUpdateDisplay(void)
{
  if (noGUI) return;

  XmUpdateDisplay(unixGUI.topLevel);

  return;
}



