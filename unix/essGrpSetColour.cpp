
#include <Xm/Xm.h>

#include "UnixGUI.h"
#include "Plot.h"
#include "FunctionsProgram.h"


extern UnixGUI unixGUI;
extern Plot    plot;
extern bool    noGUI;


void essGrpSetColour(int col, bool stdCol)
{
  if (noGUI) return;
	
  int c = col;

  if (!stdCol) 
  { 
    if (c < 0) c = 0; else if (c > 511) c = 511;
    //if (col < 0 || col > 511) prgWarning(1,"essGrpSetColour","inadmissible colour id!");
    c += plot.nNamedColours;
  }
  else
    if (col < 0 || col > plot.nNamedColours) 
      prgWarning(2,"essGrpSetColour","inadmissible colour id!");
	
  XSetForeground(XtDisplay(unixGUI.topLevel), unixGUI.gc, unixGUI.colourPixelValue[c]);

  return;
}
	
