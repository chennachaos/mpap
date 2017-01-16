
#include "Macro.h"
#include "Plot.h"
#include "RunControl.h"
#include "FunctionsEssGrp.h"
#include "DomainTree.h"
#include "Mesh.h"


extern DomainTree domain;
extern Plot       plot;
extern RunControl runCtrl;


int macro19(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "pers";
    macro.type = "plot";
    macro.what = "change perspective, 3D";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;
    macro.sensitivity[PRE]   = true;

    macro.db.selectDomain();

    macro.db.addRadioBox("*mouse","(keyboard)");
    macro.db.addToggleButton("deformed conf.",1);
    macro.db.addToggleButton("bounding box",0);
    macro.db.addToggleButton("backface swap",0);
    macro.db.addToggleButton("set des.=act.",0);
    macro.db.addToggleButton("move zoom mode",0);
    macro.db.addRadioBox("XY","YZ","ZX","*free");
    macro.db.nextButtonBox();
    macro.db.addTextField("",0,4);
    macro.db.addTextField("",0,4);
    macro.db.addTextField("",0,4);
    macro.db.addTextField("",0,4);
    macro.db.addTextField("",0,4);
    macro.db.addTextField("",0,4);
    macro.db.addTextField("",0,4);
    macro.db.addTextField("",0,4);
    macro.db.addTextField("",0,4);
    macro.db.addTextField("",0,4);
    macro.db.addTextField("",0,4);
    macro.db.addTextField("",0,4);
    macro.db.addTextField("",0,4);

    macro.db.frameRadioBox();
    macro.db.frameButtonBox();
    macro.db.setButtonColumnDim(5,3);

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int    type, id, i;
  double xmn[3], xmx[3], xmnDef[3], xmxDef[3], fact = 0.;
  bool   defm, back;

  type = roundToInt(macro.p[0]);
  id   = roundToInt(macro.p[1]) - 1;
  defm = (roundToInt(macro.p[3]) == 1);
  back = (roundToInt(macro.p[5]) == 1);

  if (domain(type,id).ndm != 3) { COUT << "'pers' works in 3D only!\n\n"; return 0; }

  RotModeEnum   oldRotMode  = plot.perspective.rotMode;
  plot.perspective.zoomMode = (ZoomModeEnum)(roundToInt(macro.p[7]));
  plot.perspective.rotMode  = (RotModeEnum) (roundToInt(macro.p[8])-1);

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//  mouse

  if (roundToInt(macro.p[2]) == 1)
  {
    if (runCtrl.mode == BATCH) 
      { COUT << "'pers' with mouse input not admissible in 'batch' mode!\n\n"; return 0; }

    COUT << "left mouse button    -> rotate\n";
    COUT << "right mouse button   -> translate\n";
    COUT << "mouse wheel          -> zomm in/out\n";
    COUT << "left + right buttons -> quit\n\n";

    plot.perspective.setDesEqAct = (roundToInt(macro.p[6]) == 1);

    domain(type,id).findMinMaxX(xmn,xmx,defm);

    if (!plot || (plot.perspective.rotMode != oldRotMode && plot.perspective.rotMode != ROTFREE)) 
    {
      plot.fit(xmn,xmx,3);
      plot.wipe();
    }

    for (i=0; i<3; i++) plot.perspective.xCntr[i] = (xmn[i] + xmx[i]) * 0.5;

    if (roundToInt(macro.p[4]) == 1)

      plot.perspective.generateBoundingBox(xmn,xmx);

    else 

      plot.perspective.changePersObjSurfPtr = mesh(domain(type,id)).surf3D;

    essGrpSetSensAllButDrawingArea(false);
   
    runCtrl.fixStatus(CHANGEPERSPECTIVE);

    plot.setColour(3);
    plot.perspective.plotChangePersObjSurf(back,defm);
  }

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// keyboard

  else
  {
    plot.perspective.prepare(macro.p,9);

    plot.dDes[0] = 1.;
    plot.dDes[1] = std::abs(macro.p[21]);

    plot.x0Des[0] = -.5;
    plot.x0Des[1] = -.5 * plot.dDes[1];

    plot.adjustToNewSize();
    essGrpWipe();
  }

//--------------------------------------------------------------------------------------------------
  return 0;  
}

