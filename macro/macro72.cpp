
#include "Macro.h"
#include "Mesh.h"
#include "DomainTree.h"


extern DomainTree domain;



int macro72(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "extr";
    macro.type = "prep";
    macro.what = "generate 3D mesh by extruding a 2D mesh";

    macro.sensitivity[PRE]   = true;
    macro.sensitivity[BATCH] = true;
    macro.sensitivity[INTER] = true;

    macro.db.selectDomain();

    macro.db.addRadioBox("x","y","*z");

    macro.db.addTextField("ndf = ",3,2);

    macro.db.addTextField("np = ",2,2);

    macro.db.nextButtonBox();
    macro.db.addTextField("p0",0,8);
    macro.db.addTextField("p1",0,8);
    macro.db.addTextField("p2",0,8);
    macro.db.addTextField("p3",0,8);
    macro.db.addTextField("h0",0,8);
    macro.db.addTextField("h1",0,8);
    macro.db.addTextField("h2",0,8);
    macro.db.addTextField("h3",0,8);

    macro.db.setButtonColumnDim(6,6,6,6);

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int typ, id, xyz, ndf, np, n;

  double sze[4], pos[4];

  typ = roundToInt(macro.p[0]);
  id  = roundToInt(macro.p[1]) - 1;

  xyz = roundToInt(macro.p[2]);

  ndf = roundToInt(macro.p[3]);
  np  = roundToInt(macro.p[4]);

  pos[0] = macro.p[5];
  pos[1] = macro.p[6];
  pos[2] = macro.p[7];
  pos[3] = macro.p[8];

  sze[0] = macro.p[ 9];
  sze[1] = macro.p[10];
  sze[2] = macro.p[11];
  sze[3] = macro.p[12];

  Domain &dom = domain(typ,id);

  cout << "   extruding " << domain.key(&dom) << " ...\n\n"; 

  domain.newDom(new Mesh);  
           
  n = domain[MESH].dom.n;

  if (!mesh(domain(MESH,n-1)).extrude(&dom,pos,sze,np,ndf))
  {
    domain.delDom(&domain(MESH,n-1));
    cout << "      mesh extrusion failed!\n\n";
    return 0;
  }
  else domain.delDom(&dom);

//--------------------------------------------------------------------------------------------------
  return 0;  
}

