
#ifndef incl_myDataIntegrateCutFEM_h
#define incl_myDataIntegrateCutFEM_h

#include "headersBasic.h"


using namespace std;
using namespace Eigen;

typedef  Vector3d  myPoint;

struct myDataIntegrateCutFEM
{
  double data[20], dvol, PENALTY, NitscheFact;
  bool isNitsche;
  int dir;
  MatrixXd  K1, K2, Kc;
  VectorXd  F1, F2;
  myPoint  param, geom, normal, specVal;

};



#endif