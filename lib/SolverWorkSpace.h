
#ifndef incl_SolverWorkSpace_h
#define incl_SolverWorkSpace_h



class SolverWorkSpace
{
  public:
	 
    SolverWorkSpace(void);

    ~SolverWorkSpace();
	  
    double *dbl;

    int  dimDbl;

    void expand(int);
    
  private:



};


#endif

