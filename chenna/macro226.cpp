
#include "Macro.h"
#include "DomainTree.h"
#include "IsogeometricFEM.h"
#include "HBSplineFEM.h"
#include "HBSplineCutFEM.h"
#include "StandardFEM.h"
#include "Plot.h"
//#include "HBScutFEMElasticity.h"


extern DomainTree domain;
extern Plot       plot;


int macro226(Macro &macro)
{
  if(!macro) 
  { 
    macro.name = "flow";
    macro.type = "vtk";
    macro.what = "contour plot for fluids";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;

    macro.db.selectDomain();


    macro.db.frameButtonBox();

    macro.db.addRadioBox("*velocity","pressure","vorticity");
    
    macro.db.frameRadioBox();

    macro.db.addRadioBox("*x-dir","y-dir","z-dir","magn");

    macro.db.addTextField("nCol ", 10);

    macro.db.addTextField("res1 ", 1);
    
    macro.db.addTextField("res2 ", 1);
    
    macro.db.addTextField("res3 ", 1);

    macro.db.frameRadioBox();

    macro.db.nextButtonBox();

    macro.db.addToggleButton("use min/max",false);
    macro.db.addTextField("min = ",0,5);
    macro.db.addTextField("max = ",1,5);

    macro.db.frameButtonBox();

    // and other stuff

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int  type, id, val1, val2, nCols, index, resln[3];
  bool  umnx;
  double  umn, umx;

  type  = roundToInt(macro.p[0]);
  id    = roundToInt(macro.p[1]) - 1;
  val1  = roundToInt(macro.p[2]) - 1;	// parameter type
  val2  = roundToInt(macro.p[3]) - 1;	// direction
  nCols = roundToInt(macro.p[4]);	// no. of colours for contour plot

  resln[0] = roundToInt(macro.p[5]) ;	// resolution to plot displacement in 1st parameter direction
  resln[1] = roundToInt(macro.p[6]) ;	// resolution to plot displacement in 2nd parameter direction
  resln[2] = roundToInt(macro.p[7]) ;	// resolution to plot displacement in 3rd parameter direction

  umnx  = (roundToInt(macro.p[8]) == 1);	// use user input min/max values
  umn   = macro.p[9];
  umx   = macro.p[10];

  if(type == 26)
    isogeometricFEM(domain(type,id)).postProcessFlow(val1, val2, nCols, umnx, umn, umx, resln);
 
  if(type == 27)
    hbsplineFEM(domain(type,id)).postProcessFlow(val1, val2, nCols, umnx, umn, umx, resln);

  if(type == 28)
    hbsplineCutFEM(domain(type,id)).postProcessFlow(val1, val2, nCols, umnx, umn, umx, resln);

  if(type == 29)
    standardFEM(domain(type,id)).postProcess(val1, val2, nCols, umnx, umn, umx, resln);

  //if(type == 30)
    //hbscutFEMElasticity(domain(type,id)).postProcessFlow(val1, val2, nCols, umnx, umn, umx, resln);

//--------------------------------------------------------------------------------------------------
  return 0;  
}

