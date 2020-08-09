

#ifndef MY_CONSTANTS_H
#define MY_CONSTANTS_H

#include <math.h>
#include <limits.h>


static double PI = acos(-1.0);
static double twoPI = 2*acos(-1.0);

//const double PI = 3.141592653589793238462643383279502884197169399375105;

#define EPSILON 1E-12
#define NEG_EPSILON  -1E-12


#define BIGD DBL_MAX
#define SMAD DBL_MIN
#define BIGI INT_MAX
#define SMAI INT_MIN

#define _NOZ +BIGD
#define _NOW -BIGD


#endif