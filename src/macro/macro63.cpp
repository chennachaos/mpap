
#include "Macro.h"
#include "DomainTree.h"
//#include "Mesh.h"
#include "RunControl.h"


extern DomainTree domain;
extern RunControl runCtrl;


int macro63(Macro &macro)
{
/*
  if (!macro) 
  { 
    macro.name = "imsh";
    macro.type = "prep";
    macro.what = "import mesh from file";

    macro.sensitivity[PRE]   = true;
    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;
    macro.sensitivity[NOPROJ]= true;

    macro.db.addList("*netgen","gmsh","mpap2/mesh","mpap2/struct2D");

    macro.db.addTextField("ndm = ",2,2);
    macro.db.addTextField("ndf = ",2,2);

    macro.db.stringTextField("file name ","mesh",40);

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int  fmt, ndm, ndf, n;

  bool isALE;

  fmt  = roundToInt(macro.p[0]);
  ndm  = roundToInt(macro.p[1]);
  ndf  = roundToInt(macro.p[2]);

  domain.newDom(new Mesh);  
           
  n = domain[MESH].dom.n;
           
  cout << "   importing MESH " << n << " ...\n\n"; 

  if (domain.ndom > 1)
    if (domain(0).ndm != ndm) 
    {
      cout << "      ndm has to match previously imported meshes.\n\n";
      cout << "      mesh import failed!\n\n"; 
      domain.delDom(&domain(MESH,n-1));
      return 0;
    }

  if (!mesh(domain(MESH,n-1)).import(fmt,macro.strg.asCharArray(),ndm,ndf))
  {
    domain.delDom(&domain(MESH,n-1));
    cout << "      mesh import failed!\n\n";
    return 0;
  }

  if (runCtrl.mode == NOPROJ) { runCtrl.newMode(PRE); runCtrl.newStatus(INTERACTIVE); }
*/
//--------------------------------------------------------------------------------------------------
  return 0;  
}

