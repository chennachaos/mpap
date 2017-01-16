#ifndef incl_ImmersedIntegrationElement_h
#define incl_ImmersedIntegrationElement_h



#include "TreeNode.h"
#include "BoundaryPoint.h"


typedef TreeNode<2>  node;
typedef BoundaryPoint<2>  Bpoint;


#ifdef _DOMAIN1D
typedef TreeNode<1>  node;
typedef BoundaryPoint<1>  Bpoint;
#endif

#ifdef _DOMAIN3D
typedef TreeNode<3>  node;
typedef BoundaryPoint<3>  Bpoint;
#endif



#endif