
#include "Macro.h"
#include "DomainTree.h"
#include "Mesh.h"


extern DomainTree domain;


static Vector<double> circ, diff;



int macro77(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "circ";
    macro.type = "wulf";
    macro.what = "tool for satisfying Kutta condition";

    macro.sensitivity[BATCH] = true;
    macro.sensitivity[INTER] = true;

    macro.db.selectDomain();
 
    macro.db.stringList("*set","chck","inc","updt");

    macro.db.addTextField("circulation / e1 ",0,6);
    macro.db.addTextField("              e2 ",0,6);
 
    macro.db.addToggleButton("print",1);

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int typ, id, i, e1, e2;

  double g1[3], g2[3];

  typ  = roundToInt(macro.p[0]);
  id   = roundToInt(macro.p[1]) - 1;

  List<DependentDoF> &uDep = mesh(domain(typ,id)).uDep;

  Element **elem = mesh(domain(typ,id)).elem;

  if (macro.strg == "chck")
  {
    e1 = roundToInt(macro.p[2]) - 1;
    e2 = roundToInt(macro.p[3]) - 1;

    elem[e1]->getGradient(1,g1);
    elem[e2]->getGradient(1,g2);

    diff.append(g1[2]-g2[2]);
 
    if (roundToInt(macro.p[4]) == 1)
    {
      COUT << "element " << e1+1 << ": grad(Phi) = (" << g1[0] << "," << g1[1] << ")\n";
      COUT << "element " << e2+1 << ": grad(Phi) = (" << g2[0] << "," << g2[1] << ")\n\n";
 
      COUT << "difference = " << diff.lastCoeff() << "\n\n";
    }

    return 0;
  }

  if (macro.strg == "set") circ.append(macro.p[2]);

  if (macro.strg == "inc")
  {
    if (circ.n < 1) prgError(1,"macro77","user error!");

    circ.append(circ.lastCoeff() + macro.p[2]);
  }
  if (macro.strg == "updt")
  {
    if (diff.n < 2) prgError(2,"macro77","user error!");

    i = circ.n - 2;

    circ.append((circ[i]*diff[i+1]-circ[i+1]*diff[i]) / (diff[i+1]-diff[i]));
  }
  for (i=0; i<uDep.n; i++)

    if (uDep[i].master.n == 1) uDep[i].ucBase = circ.lastCoeff();

  if (roundToInt(macro.p[4]) == 1) COUT << "new circulation = " << circ.lastCoeff() << "\n\n";

//--------------------------------------------------------------------------------------------------
  return 0;  
}

