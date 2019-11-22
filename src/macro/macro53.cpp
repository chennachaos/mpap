
#include "Macro.h"
//#include "Mesh.h"
#include "DomainTree.h"
#include "MathVector.h"



int macro53(Macro &macro)
{
/*
  if (!macro) 
  { 
    macro.name = "elsz";
    macro.type = "anly";
    macro.what = "calculate element sizes (prepare for 'gmsh')";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;
    macro.sensitivity[PRE]   = false;

    macro.db.selectDomain();

    macro.db.addTextField("nCircle  ",32,8);
    macro.db.addTextField("nAcross* ",4,8);
    macro.db.addTextField("maxGrad  ",0.12,8);
    macro.db.addLabel("");
    macro.db.addLabel("(*) remeshing only");

    macro.db.nextButtonBox();

    macro.db.addTextField("   hmin",0.,8);
    macro.db.addTextField("   hmax",1.0,8);
    macro.db.addTextField("errBar*",0.1,8);
    macro.db.addToggleButton("never increase el.size",false);

    macro.db.setButtonColumnDim(5,5);

    macro.db.stringTextField("select element groups* (example: 1_3_4) ","1",20);

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int  type, id, nCircle, nAcross, i, nw;

  double *x, maxGrad, hmin, hmax, errBar;

  bool keepSmallFlag;

  type    = roundToInt(macro.p[0]);
  id      = roundToInt(macro.p[1]) - 1;
  nCircle = roundToInt(abs(macro.p[2])); 
  nAcross = roundToInt(abs(macro.p[3])); 
  maxGrad = abs(macro.p[4]);        
  hmin    = abs(macro.p[5]);
  hmax    = abs(macro.p[6]);
  errBar  = macro.p[7];
  keepSmallFlag = (roundToInt(macro.p[8]) == 1);

  MyString &list = macro.strg.strip(), *word;

//  if (domain(type,id).inheritsFrom(MESH)) 
//    { COUT << "this is for Meshes and derived domains only! abort elsz.\n\n"; return 0; }

  Mesh &dom = *((Mesh*)(&(domain(type,id))));

  Vector<int> &elGrp = dom.elemGrpToBeMeshed;

  elGrp.free();

  for (i=0; i<list.length(); i++) if (list[i] == '_') list[i] = ',';

  nw = list.split(&word);

  i = 0;

  while (i < nw && word[i].toInt(elGrp.append())) i++;

  if (i < nw) { COUT << "invalid list of element groups! abort elsz.\n\n"; return 0; }

  if (elGrp.n == 0) elGrp.append(1);

  i = 0; while (i < elGrp.n) if (elGrp[i] > dom.elemGrp.n) elGrp.del(i); else i++;

  for (i=0; i<elGrp.n; i++) elGrp[i]--;

  if (!dom.geometryDiscretisedAtLeastOnce) 
    { elGrp.free(); for (i=0; i<dom.elemGrp.n; i++) elGrp.append(i); }

  dom.calcElemSizeCurr();

  if (dom.geometryDiscretisedAtLeastOnce) dom.calcElemSizeOpt(hmin,hmax,errBar,keepSmallFlag);
  else dom.getInputElemSizeOpt();

  dom.smoothElemSizeOptDist(nCircle,nAcross,maxGrad,hmin);
*/  
//--------------------------------------------------------------------------------------------------
  return 0;  
}

