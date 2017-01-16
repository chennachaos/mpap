
#ifndef incl_util_h
#define incl_util_h


#include <vector>
#include <Eigen/Dense>

#include "myConstants.h"


using Eigen::Vector3d;
using Eigen::VectorXd;
using Eigen::MatrixXd;

using std::vector;


typedef VectorXd  Point;
typedef Vector3d  point3d;
typedef Vector3d  myPoint;


// global variables for HBSpline methods

  // go generate cell edges in 2D from corners of axis-aligned bounding box 'AABB'
static  int verts[4][2] ={ {0,1}, {1,3}, {3,2}, {2,0} };
  //int verts[4][2] ={ {0,1}, {1,2}, {2,3}, {3,0} };


  // for mapping 2D faces in 3D to global 3D coordinates
static  int coord_map_HBS_2Dto3D[6][2] ={ {1,2}, {1,2}, {0,2}, {0,2}, {0,1}, {0,1} };
  //int bbox_map[3] = {};



template<typename T>
void findUnique(vector<T>& vect)
{
  sort(vect.begin(), vect.end());
  vect.erase(unique(vect.begin(), vect.end()), vect.end());
}




void GenerateCoeffMatrices(int deg, MatrixXd& SL, MatrixXd& SR);


void TensorProduct(MatrixXd& A, MatrixXd& B, MatrixXd& C);


void printMatrix(MatrixXd& AA);


void printVector(VectorXd& AA);


void printVector(vector<int>&  vec);


void printVector(vector<double>&  vec);


void printVector(double*  data, int nn);


double  HeavisideFunction(double uu, double ee);


double  DiracDelta1(double r, double alpha);


double  DiracDelta2(double r, double alpha);


double  IntegralDoubleDiracDelta1(double beta, double gamma);


void SetTimeParametersFluid(int tis, double rho, double dt, VectorXd& td);


void SetTimeParametersSolid(int tis, double rho, double dt, VectorXd& td);



void create_vector(double start, double end, double incr, vector<double>& uuu);


void map2DPointTo3DPoint(int side, myPoint& ptTemp, double val3);


// computes the factorial of a number 'nn'
// value is returned as 'double' to facilitate its direct use for 'divisions'
double factorial(unsigned int nn);


// Computes the BINOMIAL COEFFICIENT for "m" and "n"
double Bin(unsigned int m, unsigned int n);


void  getColorValue(int col, double*);


//Function to COMPARE TWO DOUBLE type values
inline  bool CompareDoubles (double A, double B) 
{
   double diff = A - B;

   return (diff < EPSILON) && (diff > NEG_EPSILON);
}


bool doubleGreater(double left, double right, bool orequal = false);


bool doubleLess(double left, double right, bool orequal = false);




inline void Idev2D(double  Idev[][4])
{
   double  r1d3 = 1.0/3.0, r2d3 = 2.0*r1d3;

   Idev[0][0] =  r2d3;     Idev[0][1] = -r1d3;	 Idev[0][2] = -r1d3;    Idev[0][3] = 0.0;
   Idev[1][0] = -r1d3;	   Idev[1][1] =  r2d3;	 Idev[1][2] = -r1d3;    Idev[1][3] = 0.0;
   Idev[2][0] = -r1d3; 	   Idev[2][1] = -r1d3;	 Idev[2][2] =  r2d3;    Idev[2][3] = 0.0;
   Idev[3][0] =  0.0;      Idev[3][1] =  0.0;    Idev[3][2] =  0.0;	Idev[3][3] = 1.0;

 return;
}

inline void Idev3D(double  Idev[][6])
{
   double  r1d3 = 1.0/3.0, r2d3 = 2.0*r1d3;

   Idev[0][0] =  r2d3;     Idev[0][1] = -r1d3;	 Idev[0][2] = -r1d3;    Idev[0][3] = 0.0;    Idev[0][4] = 0.0;    Idev[0][5] = 0.0;
   Idev[1][0] = -r1d3;	   Idev[1][1] =  r2d3;	 Idev[1][2] = -r1d3;    Idev[1][3] = 0.0;    Idev[1][4] = 0.0;    Idev[1][5] = 0.0;
   Idev[2][0] = -r1d3; 	   Idev[2][1] = -r1d3;	 Idev[2][2] =  r2d3;    Idev[2][3] = 0.0;    Idev[2][4] = 0.0;    Idev[2][5] = 0.0;
   Idev[3][0] =  0.0;      Idev[3][1] =  0.0;    Idev[3][2] =  0.0;	Idev[3][3] = 1.0;    Idev[3][4] = 0.0;    Idev[3][5] = 0.0;
   Idev[4][0] =  0.0;      Idev[4][1] =  0.0;    Idev[4][2] =  0.0;	Idev[4][3] = 0.0;    Idev[4][4] = 1.0;    Idev[4][5] = 0.0;
   Idev[5][0] =  0.0;      Idev[5][1] =  0.0;    Idev[5][2] =  0.0;	Idev[5][3] = 0.0;    Idev[5][4] = 0.0;    Idev[5][5] = 1.0;

 return;
}


double dotProductVecs(double* vec1, double* vec2, int N);


inline bool my_any_of(std::vector<int>& vecTemp, int data )
{
  bool val = false;
  for(int ii=0; ii<vecTemp.size(); ii++)
  {
    if(vecTemp[ii] == data)
    {
      val = true;
      break;
    }
  }

  return val;
}




#endif






