/*=============================================================================
        File: NurbsShapeFunctions.h
  Created by: Chennakesava Kadapa          (24 Jan 2011)
 Purpose    : Header file for the definitions of SHAPE FUNCTIONS routine 
 
 ============================================================================*/

#ifndef NurbsShapeFunctions_H
#define NurbsShapeFunctions_H


#include <iostream>
#include "NurbsUtilitiesSURFACE.h"
#include "NurbsShapeFns.h"
using namespace std;



void NurbsShapeFunctions2DAlg2(NurbsSURFACE* surf1, int ni, int nj, double u_tilde, double v_tilde, double* NN, double& J, double* dircos);


void NurbsShapeFunctions2DAlg3(NurbsSURFACE* surf1, int ni, int nj, double u, double v, double* NN, double& J);


void NurbsShapeFunctions2DAlg11(NurbsSURFACE *surf1, int ni, int nj, double u, double v, double* N, double* dN_dx, double* dN_dy, double& J);


void NurbsShapeFunctions2DAlg55(NurbsSURFACE* surf1, int ni, int nj, int confflag, double* dN_dx, double* dN_dy, double* F, double& detF);


inline void GaussPoints1D(int p, VectorArray<double>& gausspoints, VectorArray<double>& gaussweights)
{
  cerr << '\t' << " need to modify this routine " << endl;
}



double extrapolate(int ngp, int index, double* val);

double extrapolate2D(int ngp1, int ngp2, int index, double* val);

double extrapolate3D(int ngp1, int ngp2, int ngp3, int index, double* val);

bool checkClosedness(int ndim, int dir, int ngbf1, int ngbf2, int ngbf3, int* ctrlpoints);






#endif //NurbsShapeFunctions_H








