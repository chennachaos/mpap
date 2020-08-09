
#ifndef incl_FunctionsSolver_h
#define incl_FunctionsSolver_h


// functions in solver and subdirectories

extern "C"
{
  void ma38id_(int*, double*, int*);
  void myma38ad_(int*, int*, int*, int*, int*, int*, double*, int*, 
		 int*, double*, int*, int*, double*);
  void myma38bd_(int*, int*, int*, int*, int*, int*, double*, int*, 
		 int*, double*, int*, int*, double*);
  void myma38cd_(int*, int*, int*, int*, int*, double*, int*, 
		 int*, double*, double*, double*, double*, int*, int*, double*);
  
  void ma41id_(double*, int*, int*);
  void ma41ad_(int*, int*, int*, int*, int*, double*, double*, double*, double*, int*, 
	       int*, int*, double*, int*, double*, int*, int*, double*);

  void pardisoinit_(void*, int*, int*, int*, double*, int*);

  void pardiso_(void*, int*, int*, int*, int*, int*, double*,
                int*, int*, int*, int*, int*, int*, double*, double*, int*, double*);

}



#endif

