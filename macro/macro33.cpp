
#include "Macro.h"
#include "DomainTree.h"
#include "Plot.h"


extern DomainTree domain;
extern Plot       plot;


int macro33(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "1d";
    macro.type = "plot";
    macro.what = "plot one dimensional problem";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;
    macro.sensitivity[PRE]   = true;
    
    macro.db.selectDomain();

    macro.db.addList(COLOURS_BLUE);

    macro.db.addLabel       (" x-scaling");
    macro.db.addToggleButton(" current", false);
    macro.db.addTextField(" min ", 0);
    macro.db.addTextField(" max ",-1);
    
    macro.db.nextButtonBox();
    
    macro.db.addLabel       (" u-scaling");
    macro.db.addToggleButton(" current", false);
    macro.db.addTextField(" min ", 0);
    macro.db.addTextField(" max ",-1);
    
    macro.db.nextButtonBox();
    
    macro.db.addTextField   ("index ", 1,3);
    macro.db.addToggleButton("mesh",    false);
    macro.db.addToggleButton("nodes",   false);
    macro.db.addToggleButton("lines",   false);

    macro.db.addToggleButton("vertical",false);
    macro.db.addToggleButton("numbers", false);
    macro.db.addToggleButton("numbers", false);
    macro.db.addToggleButton("boundary",false);
    
    macro.db.setButtonColumnDim(4,4,4);
   
    macro.db.frameButtonBox();
    
    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int type, id, col, indx, scal;
  
  double shft, a, b, umn, umx, xmn, xmx, x1[3], x2[3],
	 xmin = 0., xmax = -1.,
	 umin = 0., umax = -1.;

  bool vert, xcur, ucur;
  
  type = roundToInt(macro.p[0]);
  id   = roundToInt(macro.p[1]) - 1;
  col  = roundToInt(macro.p[2]);
  xcur = (roundToInt(macro.p[3])==1);
  xmn  = macro.p[4];
  xmx  = macro.p[5];
  ucur = (roundToInt(macro.p[6])==1);
  umn  = macro.p[7];
  umx  = macro.p[8];
  indx = roundToInt(macro.p[9]);
  vert = (roundToInt(macro.p[13])==1);
  
  if (domain(type,id).ndm != 1) 
    {  COUT << "'plot,1d' applies to one dimensional domains only!\n\n"; return 0; }
  
  plot.setColour(col-1);

  domain(type,id).findMinMaxX(&xmin,&xmax);

  if (xmn >= xmx) { xmn = 1.1*xmin - 0.1*xmax; xmx = 1.1*xmax - 0.1*xmin; }

  if (indx > 0) 
  {
    domain(type,id).findMinMaxU(indx,umin,umax);
    if (umn >= umx) { umn = 1.2*umin - 0.2*umax; umx = 1.2*umax - 0.2*umin; }
  }
 
  plot.adjust1DSettings(xcur,xmn,xmx,ucur,umn,umx,vert);

  if (roundToInt(macro.p[10])==1) // mesh
  {
    x1[0] = xmin; x1[1] = 0.;
    x2[0] = xmax; x2[1] = 0.;
	  
    plot.line(x1,x2);
  }

  if (roundToInt(macro.p[11])==1) // nodes
  {
    if (roundToInt(macro.p[14]==1)) 
    {
    
    }
    else
    {

    }
  }

  if (roundToInt(macro.p[12])==1) // lines
  {

  }

  if (roundToInt(macro.p[15])==1) // element numbers
  {

  }

  if (roundToInt(macro.p[16])==1) // boundary lines
  {

  }

  domain(type,id).plotU1D(indx);
  
//--------------------------------------------------------------------------------------------------
  return 0;  
}

