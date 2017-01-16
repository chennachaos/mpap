/*=============================================================================
        File: NurbsUtilitiesCURVE.h
  Created by: Chennakesava Kadapa          (16 Dec 2010)
 Purpose    : Header file for the definitions of functions to modify a NURBS CURVE 
              like Knot Insertion, Degree Elevation etc...

       All functions are defined as INLINE functions 
 ============================================================================*/
#ifndef NurbsUtilitiesCURVE_H
#define NurbsUtilitiesCURVE_H


#include <iostream>
#include <math.h>
#include "NurbsCURVE.h"
using namespace std;


void CurveKnotIns(NurbsCURVE* curv1, double ubar, int r, NurbsCURVE* curv2);


void RefineCurveKnotVector(NurbsCURVE* curv1, KNOTVECTOR& X, NurbsCURVE* curv2);


void DegreeElevateCurve(NurbsCURVE* curv1, int t, NurbsCURVE* curv2);



int RepamCurveNelem(NurbsCURVE* curv1, int p, int Nel1, NurbsCURVE* curv2);


int RepamCurvenCPs(NurbsCURVE* curv1, int p, int nCP1, NurbsCURVE* curv2);



#endif  //NurbsUtilitiesCURVE_H
