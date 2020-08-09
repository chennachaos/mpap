
#include "Macro.h"
#include "DomainTree.h"
#include "HBSplineFEM.h"
#include "HBSplineCutFEM.h"
#include "StandardFEM.h"
<<<<<<< HEAD:chenna/macro201.cpp
//#include "Plot.h"
//#include "HBScutFEMElasticity.h"

extern DomainTree domain;
//extern Plot       plot;
=======

extern DomainTree domain;
>>>>>>> collabchandan:src/chenna/macro201.cpp

int macro201(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "data";
    macro.type = "chen";
    macro.what = "plot geom";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;

    macro.db.selectDomain();

    macro.db.addList(COLOURS_BLUE);
   
   
    macro.db.frameButtonBox();

    macro.db.addRadioBox("*initialgeom","finalgeom","result");
    
    macro.db.frameRadioBox();

    macro.db.frameButtonBox();

    macro.db.addRadioBox("*elements","control points.");
    
    macro.db.frameRadioBox();

    macro.db.addToggleButton("plot KNOT lines",true);

    macro.db.addTextField("res1 ", 1);
    
    macro.db.addTextField("res2 ", 1);
    
    macro.db.addTextField("res3 ", 1);

    macro.db.frameRadioBox();

    // and other stuff

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int  type, id, col, val1, resln[3];
  bool   defm = true, flag1 = true, flag2 = true, flag3 = true;
  double xmn[3], xmx[3];


  type = roundToInt(macro.p[0]);
  id   = roundToInt(macro.p[1]) - 1;
  col  = roundToInt(macro.p[2]);
  val1 = roundToInt(macro.p[3]);
  flag2 = (roundToInt(macro.p[4]) == 1);

  flag3 = (roundToInt(macro.p[5]) == 1);

  resln[0] = roundToInt(macro.p[6]) ;	// resolution to plot displacement in 1st parameter direction
  resln[1] = roundToInt(macro.p[7]) ;	// resolution to plot displacement in 2nd parameter direction
  resln[2] = roundToInt(macro.p[8]) ;	// resolution to plot displacement in 3rd parameter direction

  //plot.setColour(col-1);

/*
  if(!plot)
  {
     if(val1 == 1 || val1 == 2)
       isogeometricFEM(domain(type,id)).findMinMaxX(xmn,xmx,defm);
     else
       isogeometricFEM(domain(type,id)).findMinMaxResult(xmn,xmx,defm);

     plot.fit(xmn,xmx,domain(type,id).ndm);
  }

  isogeometricFEM(domain(type,id)).plotGeom(val1, flag2, col-1, flag3, resln);
*/

  //if(!plot)
  //{
     if(val1 == 1 || val1 == 2)
       domain(type,id).findMinMaxX(xmn,xmx,defm);
     //else
       //domain(type,id).findMinMaxResult(xmn,xmx,defm);

     //plot.fit(xmn,xmx,domain(type,id).ndm);
  //}
  
    domain(type,id).plotGeom(val1, flag2, col-1, flag3, resln);


//--------------------------------------------------------------------------------------------------
  return 0;  
}

