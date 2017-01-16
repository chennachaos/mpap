/*=============================================================================
        File: NurbsUtilitiesSOLID.h
  Created by: Chennakesava Kadapa          (09 Jan 2011)
 Purpose    : Header file for the definition of functions to modify a NURBS Surface
              like Knot Insertion, Degree Elevation etc...
              
       All functions are defined as INLINE functions 
 ============================================================================*/

#ifndef NurbsUtilitiesSOLID_H
#define NurbsUtilitiesSOLID_H


#include <iostream>
#include "NurbsSOLID.h"


using namespace std;



void RefineSolidKnotVector(NurbsSOLID* solid1, int repamID, int* tt, NurbsSOLID* solid2);



void  DegreeElevateSolid(NurbsSOLID* solid1, int* tt, NurbsSOLID* solid2);




#endif  //NurbsUtilitiesSOLID_H



