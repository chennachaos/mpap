
#ifndef incl_GUASSQUAD_h
#define incl_GUASSQUAD_h


#include "headersBasic.h"
#include "util.h"


using namespace std;
using namespace Eigen;




//template <int DIM>
class GaussQuadrature
{
    private:

        int nGP;

    public:

    //typedef Matrix<double, DIM, 1>  myPoint;

      vector<myPoint>  gausspoints;
      vector<double>  gaussweights;

      GaussQuadrature();

      ~GaussQuadrature();

      //int getDimension()
      //{ return DIM;}

      void SetGaussPoints(int nn);
      
      void reset();
};



/*
template <int DIM>
GaussQuadrature<DIM>::GaussQuadrature()
{
}



template <int DIM>
GaussQuadrature<DIM>::~GaussQuadrature()
{
}

template <int DIM>
void GaussQuadrature<DIM>::reset()
{
  gausspoints.clear();
  gaussweights.clear();
}
*/


#endif




