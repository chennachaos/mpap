
#ifndef incl_BasisFunctionsBSpline_h
#define incl_BasisFunctionsBSpline_h


//#include "headersBasic.h"



//////////////////////////////////////////////////////////////
// compute univariate B-Spline basis functions for 
// at parameter 'u'
// degree 'p'
// starting at knot 'start'
// with a knot span of 'incr'
// and store in array 'N'
//////////////////////////////////////////////////////////////

void HB_BasisFuns(int p, double start, double incr, double u, double* N);




void HB_DersBasisFuns(int p, double start, double incr, double u, int n, double** ders);




#endif
