/*=============================================================================
        File: NurbsUtilities.h
  Created by: Chennakesava Kadapa          (10 Dec 2010)
 Purpose    : Header file for the definitions of Basic Stuff
              required in creating geometry using NURBS

      Some other definitions are also given in NURBS_CPOINT_Class.h file
 ============================================================================*/
#ifndef NurbsUtilities_H
#define NurbsUtilities_H

#include <iostream>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <vector>
#include "NurbsCPOINT.h"
using namespace std;

typedef int INDEX;
typedef int DEGREE;
typedef int INTEGER;
typedef int FLAG;
typedef double REAL;

using namespace MatricesWulf;

// knot vector
typedef  VectorArray<REAL> KNOTVECTOR;


// Polygon for a Curve (Euclidean Points)
typedef  ListArray<EPOINT> POLYGON;


// Control Polygon for a Curve
typedef  ListArray<CPOINT> CPOLYGON;


// Net for a Surface (Euclidean Points)
typedef  ListArray<ListArray<EPOINT> > NET;


// Control Net for a Surface
typedef  ListArray<ListArray<CPOINT> > CNET;



// Net for a Surface (Euclidean Points)
typedef  ListArray<ListArray<ListArray<EPOINT> > > NET3D;


// Control Net for a Surface
typedef  ListArray<ListArray<ListArray<CPOINT> > > CNET3D;


void findunique(KNOTVECTOR& U, VectorArray<double>& XX);


void finduniqueInt(VectorArray<int>& U, VectorArray<int>& XX);


// creates a vector of type 'double' starting from 'start' and ending with 'end' with increments of 'incr'
void create_vector(double start, double end, double incr, VectorArray<double>& uuu);


void create_vector2(KNOTVECTOR& U, int num, KNOTVECTOR& uu1);


void GenKnotVecForRefining(KNOTVECTOR& U, int Nsub, KNOTVECTOR& X);


double CalcDist(EPOINT& P1, EPOINT& P2);


double DotProduct(EPOINT& P1, EPOINT& P2);


void CrossProduct(EPOINT& P1, EPOINT& P2, EPOINT& P3);


int binarySearch(int* sortedArray, int first, int last, int key);


int EntryExists(MatrixSparse<double>& mtx1, int& r1, int& c1);


//int findEntryIndex(MatrixSparse<double>& mtx1, int& r1, int& c1);


void SortArray(VectorArray<double>& X1);


void SortArrayInt(VectorArray<int>& X1);


double ArcTan(double x, double y);


void sub2VectorArraysInt(VectorArray<int>& X1, VectorArray<int>& X2, VectorArray<int>& XX);


void getLowerOrderKnotVector(VectorArray<double>& U, int p, VectorArray<double>& VV);


double vonMises3D(double* stre);


#endif  //NurbsUtilities_H










