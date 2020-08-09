
#ifndef incl_headersBasic_h
#define incl_headersBasic_h


#include <iostream>
#include <limits.h>
#include <float.h>
#include <vector>
#include <assert.h>
#include <math.h>
#include <cmath>
#include <stdio.h>
#include <set>
#include <fstream>
#include <string.h>
#include <iomanip>
#include <stdlib.h>
#include <algorithm>
#include <memory>

#include <numeric>
#include <iterator>
#include <functional>

using std::vector;

typedef  vector<int>  vecIntSTL;
typedef  vector<float>  vecFltSTL;
typedef  vector<double>  vecDblSTL;


enum  PhysicsType  {PHYSICS_TYPE_FLUID=0, PHYSICS_TYPE_SOLID=1};

enum  SolverType  {SOLVER_TYPE_EIGEN=0, SOLVER_TYPE_PETSC=1};

enum  {TRIA3NODE=0, QUAD4NODE, TET4NODE, HEX8NODE};




#endif


