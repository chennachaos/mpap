
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

      vector<myPoint>  gausspoints;
      vector<double>  gaussweights;

      GaussQuadrature() {}

      ~GaussQuadrature() {}

      void SetGaussPoints(int nn);
      
      void reset()
      {
        gausspoints.erase(gausspoints.begin(), gausspoints.begin()+gausspoints.size() );
        gaussweights.erase(gaussweights.begin(), gaussweights.begin()+gaussweights.size() );

        gausspoints.clear();
        gaussweights.clear();

        return;
      }
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




