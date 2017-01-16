
#include "Macro.h"
#include "Mesh.h"
#include "DomainTree.h"
#include "RunControl.h"


extern DomainTree domain;
extern RunControl runCtrl;



int macro71(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "s2dm";
    macro.type = "wulf";
    macro.what = "generate simple 2D mesh";

    macro.sensitivity[NOPROJ] = true;
    macro.sensitivity[PRE] = true;

    macro.db.addList("*rectangle", "circular arc","circular ring");

    macro.db.addList("*linear quadrilaterals","quadratic 8-noded quads",
                     "quadratic 9-noded quads",NULL);

    macro.db.addTextField("ndf = ",2,2);
    macro.db.addTextField("i1",0,8);
    macro.db.addTextField("i2",0,8);
    macro.db.addTextField("i3",0,8);
    macro.db.nextButtonBox();
    macro.db.addTextField("d1",0,6);
    macro.db.addTextField("d2",0,6);
    macro.db.addTextField("d3",0,6);
    macro.db.addTextField("d4",0,6);
    macro.db.addTextField("d5",0,6);
    macro.db.addTextField("d6",0,6);
    macro.db.addTextField("d7",0,6);
    macro.db.addTextField("d8",0,6);
    macro.db.addTextField("d9 ",0,6);
    macro.db.addTextField("d10",0,6);
    macro.db.addTextField("d11",0,6);
    macro.db.addTextField("d12",0,6);

    macro.db.setButtonColumnDim(5,4,4,4);

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int n, ndf, geom, elType, intData[5];

  double dblData[10];

  geom   = roundToInt(macro.p[0]);
  elType = roundToInt(macro.p[1]);
  ndf    = roundToInt(macro.p[2]);

  intData[0] = roundToInt(macro.p[3]);
  intData[1] = roundToInt(macro.p[4]);
  intData[2] = roundToInt(macro.p[5]);

  dblData[0] = macro.p[ 6];
  dblData[1] = macro.p[ 7];
  dblData[2] = macro.p[ 8];
  dblData[3] = macro.p[ 9];
  dblData[4] = macro.p[10];
  dblData[5] = macro.p[11];
  dblData[6] = macro.p[12];
  dblData[7] = macro.p[13];
  dblData[8] = macro.p[14];
  dblData[9] = macro.p[15];

  domain.newDom(new Mesh);  
           
  n = domain[MESH].dom.n;

  if (!mesh(domain(MESH,n-1)).generateSimple2DMesh(geom,elType,ndf,intData,dblData))
  {
    domain.delDom(&domain(MESH,n-1));
    cout << "      mesh generation failed!\n\n";
    return 0;
  }

  if (runCtrl.mode == NOPROJ) { runCtrl.newMode(PRE); runCtrl.newStatus(INTERACTIVE); }

//--------------------------------------------------------------------------------------------------
  return 0;  
}

