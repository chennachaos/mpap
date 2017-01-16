
#include <Xm/Xm.h>
#include <iostream>


#include "Definitions.h"
#include "Plot.h"
#include "UnixGUI.h"



extern Plot    plot;
extern UnixGUI unixGUI;
extern bool    noGUI;




void essGrpGetPlotAreaGeom(void)
{
  if (noGUI) return;

  plot.hPix = 0;
  
  XtVaGetValues(XtNameToWidget(unixGUI.topLevel,"*.drawingArea"), 
		XmNwidth, &(plot.wPix), XmNheight, &(plot.hPix), myNULL);

  plot.w = (double) plot.wPix;
  plot.h = (double) plot.hPix; 
  
  //COUT << "plot area width * height = " << plot.w << " pixel X " << plot.h << " pixel\n\n";
  
  return;
}

