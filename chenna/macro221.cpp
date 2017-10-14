
#include "Macro.h"
#include "DomainTree.h"


extern DomainTree domain;


int macro221(Macro &macro)
{
  if(!macro) 
  { 
    macro.name = "ctr2";
    macro.type = "vtk";
    macro.what = "continuous contour plot";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;

    macro.db.selectDomain();


    macro.db.frameButtonBox();

    macro.db.addRadioBox("*disp","total-strain","elastic-strain","plastic-strain","stress","intVar");
    
    macro.db.frameRadioBox();

    macro.db.addRadioBox("*x-dir","y-dir","z-dir","xy-dir","yz-dir","zx-dir","eqv","hpre");

    macro.db.addTextField("ind ", 1);

    macro.db.addTextField("nCol ", 10);

    macro.db.addTextField("res1 ", 5);
    
    macro.db.addTextField("res2 ", 5);
    
    macro.db.addTextField("res3 ", 5);

    macro.db.frameRadioBox();

    macro.db.nextButtonBox();

    macro.db.addToggleButton("extrap/direct",true);

    macro.db.addToggleButton("use min/max",false);
    macro.db.addTextField("min = ",0,5);
    macro.db.addTextField("max = ",1,5);

    macro.db.frameButtonBox();

    // and other stuff

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int  type, id, val1, val2, nCols, index, resln[3];
  bool  legd, umnx, flag1, defm = true;
  double xmn[3], xmx[3], umn, umx;

  xmn[0]=0; xmn[0]=0; xmx[0]=1; xmx[1]=1;

  type  = roundToInt(macro.p[0]);
  id    = roundToInt(macro.p[1]) - 1;
  val1  = roundToInt(macro.p[2]) - 1;	// parameter type
  val2  = roundToInt(macro.p[3]) - 1;	// direction
  index = roundToInt(macro.p[4]) - 1;	// index for internal variable
  nCols = roundToInt(macro.p[5]);	// no. of colours for contour plot

  resln[0] = roundToInt(macro.p[6]) ;	// resolution to plot displacement in 1st parameter direction
  resln[1] = roundToInt(macro.p[7]) ;	// resolution to plot displacement in 2nd parameter direction
  resln[2] = roundToInt(macro.p[8]) ;	// resolution to plot displacement in 3rd parameter direction

  flag1 = (roundToInt(macro.p[9]) == 1);	//extrapolate/direct flag

  umnx  = (roundToInt(macro.p[10]) == 1);	// use user input min/max values
  umn   = macro.p[11];
  umx   = macro.p[12];

  //if(!plot) 
  //{
    //domain(type,id).findMinMaxX(xmn,xmx,defm);

    //plot.fit(xmn,xmx,domain(type,id).ndm);
  //}

  //isogeometricFEM(domain(type,id)).contourplotVTK2(val1, val2, flag1, index, nCols, umnx, umn, umx, resln);

//--------------------------------------------------------------------------------------------------
  return 0;  
}

