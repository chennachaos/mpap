
#include "headersBasic.h"
#include "util.h"
#include "GaussQuadrature.h"


using namespace std;
using namespace Eigen;



GaussQuadrature::GaussQuadrature()
{
}



GaussQuadrature::~GaussQuadrature()
{
}



void GaussQuadrature::reset()
{
  gausspoints.erase(gausspoints.begin(), gausspoints.begin()+gausspoints.size() );
  gaussweights.erase(gaussweights.begin(), gaussweights.begin()+gaussweights.size() );

  gausspoints.clear();
  gaussweights.clear();
  
  //vector<myPoint>().swap(gausspoints);   // clear x reallocating
  //vector<double>().swap(gaussweights);

  return;
}





/*
template <>
void GaussQuadrature<1>::SetGaussPoints(int nn )
{
  return;
}

template <>
void GaussQuadrature<2>::SetGaussPoints(int nn )
{
  return;
}

template <>
void GaussQuadrature<3>::SetGaussPoints(int nn )
{
  return;
}
*/


