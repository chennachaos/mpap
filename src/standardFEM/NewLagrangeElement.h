
#ifndef  incl_newLagrangeElement_h
#define  incl_newLagrangeElement_h

#include "LagrangeElement.h"

#include "LagrangeElem1Dbar2Node.h"
#include "EulerBernoulliBeamElement2D.h"
#include "FrameElement2D.h"
#include "ElementGeomExactTruss2D.h"
#include "ElementGeomExactBeam2D.h"
#include "LagrangeElem2DStructSolidTria3Node.h"
#include "LagrangeElem2DStructSolidQuad4Node.h"
#include "LagrangeElem2DStructSolidMixed.h"
#include "LagrangeElem2DBbarFbar.h"
#include "EulerBernoulliBeamElement3D.h"
#include "LagrangeElem3DStructSolidTet4Node.h"
#include "LagrangeElem3DStructSolidHex8Node.h"
#include "LagrangeElem3DStructSolidMixed.h"
#include "LagrangeElem3DBbarFbar.h"
#include "LagrangeElem2DStructSolidTria3NodeStab.h"
#include "LagrangeElem2DStructSolidQuad4NodeStab.h"
#include "LagrangeElem3DStructSolidTet4NodeStab.h"
#include "LagrangeElem3DStructSolidHex8NodeStab.h"

#include "MindlinPlateElement.h"
#include "KirchhoffPlateElement.h"
#include "LagrangeElem3DShellQuad4Node.h"

#include "ContactElement2D1nodedContactAlongXaxis.h"
#include "ContactElement2D1nodedContactAlongYaxis.h"

#include "ContactElement3D1nodedContactAlongXaxis.h"
#include "ContactElement3D1nodedContactAlongYaxis.h"
#include "ContactElement3D1nodedContactAlongZaxis.h"


inline  LagrangeElement*  NewLagrangeElement(int type)
{
  switch (type)
  {
    case  0: return (LagrangeElement*) new LagrangeElem1Dbar2Node; break;

    case  1: return (LagrangeElement*) new EulerBernoulliBeamElement2D; break;

    case  2: return (LagrangeElement*) new FrameElement2D; break;

    case  3: return (LagrangeElement*) new ElementGeomExactTruss2D; break;

    case  4: return (LagrangeElement*) new ElementGeomExactBeam2D; break;

    case  7: return (LagrangeElement*) new LagrangeElem2DStructSolidTria3Node; break;

    case  8: return (LagrangeElement*) new LagrangeElem2DStructSolidQuad4Node; break;

    case  9: return (LagrangeElement*) new LagrangeElem2DStructSolidMixed; break;

    case 10: return (LagrangeElement*) new LagrangeElem2DBbarFbar; break;

    case 15: return (LagrangeElement*) new EulerBernoulliBeamElement3D; break;

    case 18: return (LagrangeElement*) new LagrangeElem3DStructSolidTet4Node; break;

    case 19: return (LagrangeElement*) new LagrangeElem3DStructSolidHex8Node; break;

    case 20: return (LagrangeElement*) new LagrangeElem3DStructSolidMixed; break;

    case 21: return (LagrangeElement*) new LagrangeElem3DBbarFbar; break;

    case 26: return (LagrangeElement*) new LagrangeElem2DStructSolidTria3NodeStab; break;

    case 27: return (LagrangeElement*) new LagrangeElem2DStructSolidQuad4NodeStab; break;

    case 28: return (LagrangeElement*) new LagrangeElem3DStructSolidTet4NodeStab; break;

    case 29: return (LagrangeElement*) new LagrangeElem3DStructSolidHex8NodeStab; break;

    case 30: return (LagrangeElement*) new MindlinPlateElement; break;

    case 31: return (LagrangeElement*) new KirchhoffPlateElement; break;

    case 32: return (LagrangeElement*) new LagrangeElem3DShellQuad4Node; break;

    case 33: return (LagrangeElement*) new ContactElement2D1nodedContactAlongXaxis; break;

    case 34: return (LagrangeElement*) new ContactElement2D1nodedContactAlongYaxis; break;

    case 35: return (LagrangeElement*) new ContactElement3D1nodedContactAlongXaxis; break;

    case 36: return (LagrangeElement*) new ContactElement3D1nodedContactAlongYaxis; break;

    case 37: return (LagrangeElement*) new ContactElement3D1nodedContactAlongZaxis; break;

    default: prgError(1,"ImmersedFlexibleSolid::newElement","unknown element type name!"); return NULL;
  }
}


#endif

