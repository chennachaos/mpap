/*=============================================================================
        File: NurbsUtilitiesSURFACE.h
  Created by: Chennakesava Kadapa          (09 Jan 2011)
 Purpose    : Header file for the definition of functions to modify a NURBS Surface
              like Knot Insertion, Degree Elevation etc...
              
       All functions are defined as INLINE functions 
 ============================================================================*/

#ifndef NurbsUtilitiesSURFACE_H
#define NurbsUtilitiesSURFACE_H


#include <iostream>
#include "NurbsSURFACE.h"


using namespace std;




void RefineSurfKnotVector1D(NurbsSURFACE* surf1, KNOTVECTOR& X, int dir, NurbsSURFACE* surf2);


void DegreeElevateSurf1D(NurbsSURFACE* surf1, int t, int dir, NurbsSURFACE* surf2);


void FindKnotVector(NurbsSURFACE* surf1, int dir, int Nsub, KNOTVECTOR& XX);


int RepamSurf2DNelem(NurbsSURFACE* surf0, int p, int q, int Nel1, int Nel2, NurbsSURFACE* surf1);


int RepamSurf2DnCPs(NurbsSURFACE* surf0, int p, int q, int nCP1, int nCP2, NurbsSURFACE* surf1);





/*
  Refine the Surface knot vector in both the directions
  based on Algorithm A5.5 on pg#167 of the Nurbs Book
  
  X the new knots to insert in the knot vector
*/
/*
inline void RefineSurfKnotVector2D(NurbsSURFACE* surf1, KNOTVECTOR& XU, KNOTVECTOR& XV, NurbsSURFACE* surf2)
{
    //refine the knot vector in U direction
    NurbsSURFACE* surf3 = surf1;
       
    RefineSurfKnotVector1D(surf1, XU, 1, surf3);

    //refine the knot vector in V direction
       
    RefineSurfKnotVector1D(surf3, XV, 2, surf2);
}
*/


/*
inline void DegreeElevateSurf2D(NurbsSURFACE* surf1, int t1, int t2, NurbsSURFACE* surf2)
{
    //degree elevate the surface in U direction
    NurbsSURFACE* surf_temp = surf1;

    DegreeElevateSurf1D(surf1, t1, 1, surf_temp);
       
    //degree elevate the surface in V direction
     
    DegreeElevateSurf1D(surf_temp, t2, 2, surf2);
}
*/






#endif  //NurbsUtilitiesSURFACE_H



