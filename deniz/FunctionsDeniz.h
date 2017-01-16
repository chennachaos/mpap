
#ifndef incl_FunctionsDeniz_h
#define incl_FunctionsDeniz_h
#include "Element.h"
#include "MathMatrix.h"


// functions in deniz

    int  micro4to1 (int , int , int , int, int);

    int  micro2to1 (int , int , int );

   int integrateQuadraticBoundary(double*, double*, double*, int, int);

   int integrateLinearBoundary(double*, double*, double*, int, int);



extern "C"
{

 //Material

 void multiscalematerial2d_(double*,double*,double*,double*,
		             int*,int*,int*,
		             Element*);
  
  void multiscalematerial3d_(double*,double*,double*,double*,
		             int*,int*,int*,
		             Element*);


  void hyplasmaterial_(double*,double*,double*,
					 double*,double*,double*,double*,
		             int*,int*,int*,int*,int*,
		             Element*);


  void hyplas_(double*,double*,double*,
					double*,double*,double*,double*,double*,
					int*,int*,int*,int*,int*,int*,int*);


  //Support
  void compxigpt_(double*, double*, int*);
  
  void compshpt_ (double*, double*, double*, double*, int*);
    void compshpbt_ (double*, double*, double*, double*, int*);
	  void compshpbq_ (double*, double*, double*, double*, int*);
	  void inverse_matrix_(double*,double*,int*,int*);
	 int decomplr_matrix_(double*,int*,int*,int*);

}




#endif

