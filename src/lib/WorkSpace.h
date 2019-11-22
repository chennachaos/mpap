
#ifndef incl_WorkSpace_h
#define incl_WorkSpace_h



class WorkSpace
{
  public:
	 
    WorkSpace(void);

    ~WorkSpace();
	  
    double *dblSpc;
    int    *intSpc;

    int  dimDbl, dimInt;

    void expandDbl(int);
    void expandInt(int);
    
  private:



};


#endif

