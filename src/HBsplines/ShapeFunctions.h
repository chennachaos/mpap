
#ifndef incl_ShapeFunctions_CLASS_h
#define incl_ShapeFunctions_CLASS_h

#include <Eigen/Dense>

using Eigen::VectorXd;


class ShapeFunctions
{
  public:

    int  DIM;
    
    VectorXd  N, dN_dx, d2N_dx2, dN_dy, d2N_dy2, dN_dz, d2N_dz2;

    ShapeFunctions() {}

    void initialise(int p);
    
    void initialise(int p, int q);
    
    void initialise(int p, int q, int r);

    ~ShapeFunctions() {}
    
    void setDimension(int d)
    {
      DIM = d;
    }

    void Compute();

    void zero();
    
};



#endif

