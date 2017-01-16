
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
  gausspoints.clear();
  gaussweights.clear();
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


