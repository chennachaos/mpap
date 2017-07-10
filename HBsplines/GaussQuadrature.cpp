
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

  return;
}



