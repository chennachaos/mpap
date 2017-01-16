
#ifndef incl_FunctionsMaterial_h
#define incl_FunctionsMaterial_h


#include "Element.h"


// functions in material


void matDiffTangentTest(double, char*, int, int, bool);


extern "C"
{
  int  matdim_(int*);
	
  void matlib1d_(double*,double*,double*,double*,double*,double*,double*,double*,
		 int*,int*,int*,int*,int*,int*,
		 int*      gp = NULL,
		 Element* elm = NULL);
  
  void matlib2d_(double*,double*,double*,double*,double*,double*,double*,double*,
		 int*,int*,int*,int*,int*,int*,
 		 int*      gp = NULL,
		 Element* elm = NULL);
 
  void matlib3d_(double*,double*,double*,double*,double*,double*,double*,
		 int*,int*,int*,int*,int*,
		 int*      gp = NULL,
		 Element* elm = NULL);

  void calcb_(double*, double*);
  void calcc_(double*, double*);

  void taub2sigccpart1_(double*, double*, double*, double*, double*, double*);
  void taub2sigccpart2_(double*, double*);

}




#endif

