
#include "Macro.h"
//#include "Plot.h"
#include "RunControl.h"
#include "FunctionsEssGrp.h"


////extern Plot       plot;
extern RunControl runCtrl;


int macro18(Macro &macro)
{
/*
  if (!macro) 
  { 
    macro.name = "zoom";
    macro.type = "plot";
    macro.what = "2D zoom in / out";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;
    macro.sensitivity[PRE]   = true;

    macro.db.addRadioBox("*mouse","keyboard");
    macro.db.addTextField("x0 = ");
    macro.db.addTextField("y0 = ");
    macro.db.addTextField("w = ");
    macro.db.addTextField("h = ");
    macro.db.setButtonColumnDim(2);
    macro.db.frameRadioBox();
    macro.db.frameButtonBox();
    
    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  if (plot.dim != 2) { COUT << "'zoom' works in 2D only!\n\n"; return 0; }

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//  mouse

  if (roundToInt(macro.p[0]) == 1)
  {	  
    if (runCtrl.mode == BATCH) 
      { COUT << "'zoom' with mouse input not admissible in 'batch' mode!\n\n"; return 0; }
	 
    essGrpSetSensAllButDrawingArea(false);
   
    runCtrl.fixStatus(SELECTZOOMBOX);
  }

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// keyboard

  else
  {
    if (macro.p[3] < 1.e-20 || macro.p[4] < 1.e-20) 
    {
      COUT << "width and/or height < 0 ?\n\n";
      return 0;
    }
	  
    plot.x0Des[0] = macro.p[1];
    plot.x0Des[1] = macro.p[2];
	    
    plot.dDes[0] = macro.p[3];
    plot.dDes[1] = macro.p[4];
   
    plot.adjustToNewSize(); 
  }
*/
//--------------------------------------------------------------------------------------------------
  return 0;  
}

