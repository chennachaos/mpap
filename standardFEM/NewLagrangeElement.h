
#ifndef  incl_newLagrangeElement_h
#define  incl_newLagrangeElement_h

#include "LagrangeElement.h"

#include "LagrangeElem1Dbar2Node.h"
#include "EulerBernoulliBeamElement2D.h"
#include "FrameElement2D.h"
#include "ElementGeomExactTruss2D.h"
#include "ElementGeomExactBeam2D.h"
#include "LagrangeElem2DPoissonTria3Node.h"
#include "LagrangeElem2DPoissonQuad4Node.h"
#include "LagrangeElem2DStructSolidTria3Node.h"
#include "LagrangeElem2DStructSolidQuad4Node.h"
#include "LagrangeElem2DStructSolidMixed.h"
#include "LagrangeElem2DBbarFbar.h"
#include "LagrangeElem2DStokesTria3Node.h"
#include "LagrangeElem2DStokesQuad4Node.h"
#include "LagrangeElem2DNavierStokesTria3Node.h"
#include "LagrangeElem2DNavierStokesQuad4Node.h"
#include "EulerBernoulliBeamElement3D.h"
#include "LagrangeElem3DPoissonTet4Node.h"
#include "LagrangeElem3DPoissonHex8Node.h"
#include "LagrangeElem3DStructSolidTet4Node.h"
#include "LagrangeElem3DStructSolidHex8Node.h"
#include "LagrangeElem3DStructSolidMixed.h"
#include "LagrangeElem3DBbarFbar.h"
#include "LagrangeElem3DStokesTet4Node.h"
#include "LagrangeElem3DStokesHex8Node.h"
#include "LagrangeElem3DNavierStokesTet4Node.h"
#include "LagrangeElem3DNavierStokesHex8Node.h"
#include "LagrangeElem2DStructSolidTria3NodeStab.h"
#include "LagrangeElem2DStructSolidQuad4NodeStab.h"
#include "LagrangeElem3DStructSolidTet4NodeStab.h"
#include "LagrangeElem3DStructSolidHex8NodeStab.h"

#include "MindlinPlateElement.h"
#include "KirchhoffPlateElement.h"
#include "LagrangeElem3DShellQuad4Node.h"


inline  LagrangeElement*  NewLagrangeElement(int type)
{
  switch (type)
  {
    case  0: return (LagrangeElement*) new LagrangeElem1Dbar2Node; break;

    case  1: return (LagrangeElement*) new EulerBernoulliBeamElement2D; break;

    case  2: return (LagrangeElement*) new FrameElement2D; break;

    case  3: return (LagrangeElement*) new ElementGeomExactTruss2D; break;

    case  4: return (LagrangeElement*) new ElementGeomExactBeam2D; break;

    case  5: return (LagrangeElement*) new LagrangeElem2DPoissonTria3Node; break;

    case  6: return (LagrangeElement*) new LagrangeElem2DPoissonQuad4Node; break;

    case  7: return (LagrangeElement*) new LagrangeElem2DStructSolidTria3Node; break;

    case  8: return (LagrangeElement*) new LagrangeElem2DStructSolidQuad4Node; break;

    case  9: return (LagrangeElement*) new LagrangeElem2DStructSolidMixed; break;

    case 10: return (LagrangeElement*) new LagrangeElem2DBbarFbar; break;

    case 11: return (LagrangeElement*) new LagrangeElem2DStokesTria3Node; break;
    
    case 12: return (LagrangeElement*) new LagrangeElem2DStokesQuad4Node; break;

    case 13: return (LagrangeElement*) new LagrangeElem2DNavierStokesTria3Node; break;
    
    case 14: return (LagrangeElement*) new LagrangeElem2DNavierStokesQuad4Node; break;

    case 15: return (LagrangeElement*) new EulerBernoulliBeamElement3D; break;

    case 16: return (LagrangeElement*) new LagrangeElem3DPoissonTet4Node; break;

    case 17: return (LagrangeElement*) new LagrangeElem3DPoissonHex8Node; break;

    case 18: return (LagrangeElement*) new LagrangeElem3DStructSolidTet4Node; break;

    case 19: return (LagrangeElement*) new LagrangeElem3DStructSolidHex8Node; break;

    case 20: return (LagrangeElement*) new LagrangeElem3DStructSolidMixed; break;

    case 21: return (LagrangeElement*) new LagrangeElem3DBbarFbar; break;

    case 22: return (LagrangeElement*) new LagrangeElem3DStokesTet4Node; break;

    case 23: return (LagrangeElement*) new LagrangeElem3DStokesHex8Node; break;

    case 24: return (LagrangeElement*) new LagrangeElem3DNavierStokesTet4Node; break;

    case 25: return (LagrangeElement*) new LagrangeElem3DNavierStokesHex8Node; break;

    case 26: return (LagrangeElement*) new LagrangeElem2DStructSolidTria3NodeStab; break;

    case 27: return (LagrangeElement*) new LagrangeElem2DStructSolidQuad4NodeStab; break;

    case 28: return (LagrangeElement*) new LagrangeElem3DStructSolidTet4NodeStab; break;

    case 29: return (LagrangeElement*) new LagrangeElem3DStructSolidHex8NodeStab; break;

    case 30: return (LagrangeElement*) new MindlinPlateElement; break;

    case 31: return (LagrangeElement*) new KirchhoffPlateElement; break;

    case 32: return (LagrangeElement*) new LagrangeElem3DShellQuad4Node; break;

    default: prgError(1,"ImmersedFlexibleSolid::newElement","unknown element type name!"); return NULL;
  }
}


#endif

