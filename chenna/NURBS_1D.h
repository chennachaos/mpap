/*=============================================================================
        File : NURBS_1D.h
  Created by : Chennakesava Kadapa          (16 Dec 2010)
  Purpose    : Header file for the definition of functions related to 
		 Univariate B-Splines.

       All functions are defined as INLINE functions 
 ============================================================================*/
#ifndef NURBS_1D_H
#define NURBS_1D_H

#include <iostream>
#include "NurbsUtilities.h"
using namespace std;



bool IsKnot(double* U, int Un, double u);


int FindSpan(double* U, int Un, int deg, double u);


int FindMult(double* U, int Un, DEGREE deg, double ubar);


double OneBasisFun_recurs(KNOTVECTOR& U, DEGREE deg, int i, double u);


double DersOneBasisFun_recurs(KNOTVECTOR& U, DEGREE deg, int i, double u, int der_order);


void BasisFuns(double* U, int Un, DEGREE p, double u, double* N);


void BasisFuns2D(double* U, int Un, int p, double* V, int Vn, int q, double u, double v, double* NN);


void DersBasisFuns(double* U, int Un, DEGREE p, double u, int n, double** ders);


double OneBasisAlg1(double* U, int Un, DEGREE p, int i, double u);

//void DersBasisFuns(KNOTVECTOR& U, DEGREE p, double u, int n, ListArray<VectorArray<double> >& ders);

//int FindSpan(KNOTVECTOR& U, DEGREE deg, double u);


#endif  //NURBS_1D_H
