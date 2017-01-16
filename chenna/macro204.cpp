
#include "Macro.h"
#include "DomainTree.h"
#include "IsogeometricFEM.h"
#include "Plot.h"


extern DomainTree domain;
extern Plot       plot;


int macro204(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "resu";
    macro.type = "chen";
    macro.what = "continuous contour plot";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;

    macro.db.selectDomain();


    macro.db.frameButtonBox();

    macro.db.addRadioBox("*disp","total-strain","elastic-strain","plastic-strain","stress","intVar");
    
    macro.db.frameRadioBox();

    macro.db.addRadioBox("*x-dir","y-dir","z-dir","xy-dir","eqv","hpre");

    macro.db.addTextField("index = ",1);

    macro.db.addTextField("nColors = ",10);

//    macro.db.frameRadioBox();

    macro.db.addToggleButton("extrapolate/direct",true);

    macro.db.nextButtonBox();

    macro.db.addToggleButton("show legend",true);
    macro.db.addToggleButton("use min/max",false);
    macro.db.addTextField("min = ",0,5);
    macro.db.addTextField("max = ",1,5);

    
    macro.db.frameButtonBox();

    // and other stuff

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int  type, id, val1, val2, nCols, index;
  bool  legd, umnx, flag1, defm = true;
  double xmn[3], xmx[3], umn, umx;

  xmn[0]=0; xmn[0]=0; xmx[0]=1; xmx[1]=1;

  type  = roundToInt(macro.p[0]);
  id    = roundToInt(macro.p[1]) - 1;
  val1  = roundToInt(macro.p[2]) - 1;	// parameter type
  val2  = roundToInt(macro.p[3]) - 1;	// direction
  index = roundToInt(macro.p[4]) - 1;	// index for internal variable
  nCols = roundToInt(macro.p[5]) - 1;	// no. of colours for contour plot

  flag1 = (roundToInt(macro.p[6]) == 1);	//extrapolate/direct flag

  legd = (roundToInt(macro.p[7]) == 1);	// show contour legend
  umnx  = (roundToInt(macro.p[8]) == 1);	// use user input min/max values
  umn   = macro.p[9];
  umx   = macro.p[10];


  if(!plot) 
  {
     domain(type,id).findMinMaxX(xmn,xmx,defm);

     plot.fit(xmn,xmx,domain(type,id).ndm);
  }

  isogeometricFEM(domain(type,id)).contourplot(val1, val2, flag1, index, nCols, umnx, umn, umx, legd);


  
  


//--------------------------------------------------------------------------------------------------
  return 0;  
}

