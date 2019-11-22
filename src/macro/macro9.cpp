
#include "Macro.h"
#include "DomainTree.h"
//#include "Plot.h"


////extern Plot plot;


int macro9(Macro &macro)
{
/*
  if (!macro) 
  { 
    macro.name = "cont";
    macro.type = "plot";
    macro.what = "colour contour plot";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;
    macro.sensitivity[PRE] = true;    // to be removed!!!
    
    macro.db.selectDomain();
   
    macro.db.addRadioBox("*degree of freedom",
		         "initial coordinates",
		         "current coordinates",
			 "last projection");
    
    macro.db.addTextField("index   = ",1,5);
    macro.db.addTextField("colours = ",11,5);

    macro.db.addToggleButton("show legend",false);
    
    macro.db.frameRadioBox();

    macro.db.nextButtonBox();

    macro.db.addToggleButton("use min/max",false);
    macro.db.addTextField("min = ",0,5);
    macro.db.addTextField("max = ",1,5);
    
    macro.db.frameButtonBox();

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

    int type, id, var, i, nCol;

    double umn, umx, xmn[3], xmx[3];

    bool legd, umnx;
    
    type = roundToInt(macro.p[0]);
    id   = roundToInt(macro.p[1]) - 1;
    var  = roundToInt(macro.p[2]);
    i    = roundToInt(macro.p[3]);
    nCol = roundToInt(macro.p[4]);
    legd = (roundToInt(macro.p[5]) == 1);
    umnx = (roundToInt(macro.p[6]) == 1);
    umn  = macro.p[7];
    umx  = macro.p[8];

    if (!plot) 
    {
       domain(type,id).findMinMaxX(xmn,xmx,true);

       plot.fit(xmn,xmx,domain(type,id).ndm);
    }
    
    domain(type,id).contourPlot(var,i,nCol,umnx,umn,umx,legd);
*/    
//--------------------------------------------------------------------------------------------------
  return 0;  
}

