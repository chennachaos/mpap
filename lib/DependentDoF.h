
#ifndef incl_DependentDoF_h
#define incl_DependentDoF_h

#include "MathVector.h"
#include "List.h"


class DependentDoF: public ListItem
{ 
  public:

    DependentDoF(void);
    virtual ~DependentDoF();
    
    int    nd, dof;
    
    double ucBase, duc, uc;

    int    tmFct;
    
    VectorArray<int>    master, masterNd, masterDoF;

    VectorArray<double> alpha, beta;

    void init(void);
    
    void timeUpdate(void);
   
    void update(double*, int&);
    
    void updateWithIncrement(double*, int&);

    void printInfo(void);
    
  private:

};

#endif




