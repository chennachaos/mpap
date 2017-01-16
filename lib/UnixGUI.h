
#ifndef incl_UnixGUI_h
#define incl_UnixGUI_h


#include <Xm/Xm.h>


#define myNULL (void*)0



struct UnixGUI
{
  UnixGUI(void) { XPoints = new XPoint[100]; 
	          nXPoints = 100; 
		  return; }

  ~UnixGUI()    { delete [] XPoints; return; }
	
  XtAppContext app;

  Pixmap pixmap;

  Widget topLevel, drawingArea;

  GC gc;

  XPoint *XPoints;
  
  int nXPoints;
  
  int width, height, widthMM, heightMM;

  int colourPixelValue[1000];

  XFontStruct *fontStruct;
};


#endif

