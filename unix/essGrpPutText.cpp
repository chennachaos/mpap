
#include <Xm/Xm.h>
#include <iostream>


#include "MathBasic.h"
#include "UnixGUI.h"


extern UnixGUI unixGUI;
extern bool    noGUI;


void essGrpPutText(int *x0, char* strg, int posFlg, bool BGFlag)
{
  if (noGUI) return;
	
  int l = strlen(strg), w = XTextWidth(unixGUI.fontStruct,strg,l), 
      x, y, h = unixGUI.fontStruct->ascent, wd2 = intDiv(w,2), hd2 = intDiv(h,2);

  switch (posFlg)
  {
    case  2: x = x0[0]-wd2; y = x0[1]    ; break;  // bottom centre
    case  3: x = x0[0]-  w; y = x0[1]    ; break;  // bottom right
    case  4: x = x0[0]    ; y = x0[1]+hd2; break;  // centre left
    case  5: x = x0[0]-wd2; y = x0[1]+hd2; break;  // centre centre
    case  6: x = x0[0]-  w; y = x0[1]+hd2; break;  // centre right
    case  7: x = x0[0]    ; y = x0[1]+  h; break;  // top    left
    case  8: x = x0[0]-wd2; y = x0[1]+  h; break;  // top    centre
    case  9: x = x0[0]-  w; y = x0[1]+  h; break;  // top    right
    default: x = x0[0]    ; y = x0[1]    ; break;  // bottom left
  }

  if (!BGFlag) 
    XDrawString(XtDisplay(unixGUI.topLevel), unixGUI.pixmap, unixGUI.gc, x, y, strg, l);	  
  else 
    XDrawImageString(XtDisplay(unixGUI.topLevel), unixGUI.pixmap, unixGUI.gc, x, y, strg, l);

//  XDrawString(XtDisplay(unixGUI.topLevel), XtWindow(unixGUI.drawingArea), unixGUI.gc,  
//	      x, y, strg, l);
	
  return;
}




