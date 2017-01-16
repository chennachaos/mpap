
#ifndef incl_FunctionsElement_h
#define incl_FunctionsElement_h


#include "Element.h"


// functions in elements

extern "C" 
{ 

  int element1dsolidale_(double*,double*,double*,double*,double*,double*);

  int element2dsolid_(double*,double*,double*,double*,double*,double*,double*,
		              double*,double*,
			      int*,int*,int*,int*,int*,int*,int*,int*,int*,
			      Element* elm = NULL);

  int element2dfbarsolid_(double*,double*,double*,double*,double*,double*,double*,
		                double*,double*,
			        int*,int*,int*,int*,int*,int*,int*,int*,int*,
			        Element* elm = NULL);

  int element3dsolid_(double*,double*,double*,double*,double*,double*,double*,
		               double*,double*,
			       int*,int*,int*,int*,int*,int*,int*,int*,int*,
			       Element* elm = NULL);

  int element3dfbarsolid_(double*,double*,double*,double*,double*,double*,double*,
		               double*,double*,
			       int*,int*,int*,int*,int*,int*,int*,int*,int*,
			       Element* elm = NULL);
  
  int element2dgeomexsmallstrainbeam_(double*,double*,double*,double*,double*,double*,
                                      double*,double*,double*,
				      int*,int*,int*,int*,int*,int*,int*,int*,
			              Element* elm = NULL);

  int element2dkirchhoffbeam_(double*,double*,double*,double*,double*,double*,double*,
		              double*,double*,
			      int*,int*,int*,int*,int*,int*,int*,int*,
			      Element* elm = NULL);

  int element2dtruss_(double*,double*,double*,double*,double*,double*,double*,double*,double*,
		      double*,int*,int*,int*,int*,int*,int*,int*,int*,
                      Element* elm = NULL);


  int element2dstabincompfluid_(double*,double*,double*,double*,double*,double*,double*,
		                int*,int*,int*,int*,int*);

  int element2dstabincomphighrefluid_(double*,double*,double*,double*,double*,double*,double*,
		                      int*,int*,int*,int*,int*);

  int element2dstabcompfluid_(double*,double*,double*,double*,double*,double*,double*,
		              int*,int*,int*,int*,int*);


  int element3dstabincompfluid_(double*,double*,double*,double*,double*,double*,double*,
		                int*,int*,int*,int*,int*);

  int ale2d3nodedtrianglecellcentroid_    (double*, double*, double*, double*, int*, int*, int*);
  int ale2d3nodedtriangleaspectratio_     (double*, double*, double*, double*, int*, int*, int*);
  int ale2d3nodedtrianglelinearelasticity_(double*, double*, double*, double*, int*, int*, int*);
  int ale2d3nodedtrianglehyperelasticity_ (double*, double*, double*, double*, int*, int*, int*);
  int ale2d4nodedquadrilateralaspectratio_(double*, double*, double*, double*, int*, int*, int*);

  int ale3d4nodedtetrahedronaspectratio_  (double*, double*, double*, double*, int*, int*, int*);
  int ale3d4nodedtetrahedronaspectratiotest_(double*, double*, double*, double*, int*, int*, int*);
  int ale3d4nodedtetrahedronhyperelasticity_(double*, double*, double*, double*, int*, int*);

  void compxigp2d_(double*, double*, int*);

  void compshp2d_ (double*, double*, double*, double*, int*);

  void compxigp3d_(double*, double*, int*);
  
  void compshp3d_ (double*, double*, double*, double*, int*);
  
  void trishp_(double*, double*, int*, int*, double*, double*);

}



#endif

