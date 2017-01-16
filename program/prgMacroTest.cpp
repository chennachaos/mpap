
#include "MathVector.h"
#include "Plot.h"
#include "DomainTree.h"
#include "FunctionsEssGrp.h"
#include "MacroQueue.h"



extern DomainTree domain;
extern Plot       plot;
extern MacroQueue macroQueue;



void prgMacroTest(double xPos, double yPos)
{
  int type = roundToInt(macroQueue.p[0]),
      id   = roundToInt(macroQueue.p[1]) - 1,
      isw  = roundToInt(macroQueue.p[2]);

  switch (isw)
  {
    case  1: if (domain(type,id).ndm != 2) return;

             double x[2];

             x[0] =  xPos            / plot.wFactAct + plot.x0Act[0];
             x[1] = (plot.hPix-yPos) / plot.hFactAct + plot.x0Act[1];

             plot.setColour(3);
             plot.point(x);

             essGrpCopyPixmap();

             domain(type,id).transferNodalDataTest(x);

             break;

    default: COUT << "what shall we do with the drunken sailor?\n\n";
  }

  return;
}

