
#ifndef incl_GeomDataHBSplines_h
#define incl_GeomDataHBSplines_h

//#include "headersBoost.h"

#include "ImmersedSolid.h"
#include "ShapeFunctions.h"

#include "GaussQuadrature.h"

//#include "myPolyList.h"

using namespace myGeom;


class  Function;
class  DistanceFunction;


class  GeomDataHBSplines
{
  private:
    int  DIM, degree[3], nelem[3], nlbf, ngbf[3], nGP[3], totalNGP, ndof, tis;
    double  origin[3], dX[3], gridLEN[3], Jacobian[3], Jfull, rho;

    //int ROWS;
    //double  **ders1, **ders2, **ders3;

  public:

    int firstIter;

    vector<bool>  domainInclYesNo, domainFixedYesNo;

    vector<double>  gausspoints1, gaussweights1, gausspoints2, gaussweights2, gausspoints3, gaussweights3;
    //vector<double>  boundaryGPs1, boundaryGWs1, boundaryGPs2, boundaryGWs2, boundaryGPs3, boundaryGWs3;
    vector<double>  gpsLeft, gpsRight, gwsLeft, gwsRight, FluidProps, cutFEMparams;

    vector<myPoint>  gausspoints, gausspoints2D, boundaryNormals;
    vector<double>  gaussweights, gaussweights2D;
    
    vector<vector<double> >  boundaryJacobians;
    
    vector<GaussQuadrature>  boundaryQuadrature2D;

    vector<GaussQuadrature>  boundaryQuadrature3D;

    VectorXd  bodyForce, bodyForceCur, bodyForcePrev, bodyForcePrev2, bodyForcePred, bodyForcePredPrev, bodyForcePredPrev2;

    MatrixXd  coeffLeft, coeffRight, coeffSW, coeffNW, coeffSE, coeffNE;
    MatrixXd  coeffSW_Front, coeffNW_Front, coeffSE_Front, coeffNE_Front;
    MatrixXd  coeffSW_Back, coeffNW_Back, coeffSE_Back, coeffNE_Back;

    vector<int>  node_map_new_to_old;
    vector<int>  node_map_old_to_new;


    vector<vector<ShapeFunctions> >  shpfns;

    //vector<vector<Function> >  DirichletBoundaryFunctions;
    Function*  analyDBC;

    //vector<bgpolygon_type>  polyImm;

    //myPolyList  polyImm;

    vector<DistanceFunction*>  distFuncs;
    
    vector<ImmersedSolid*>  immSolidPtrs;

    GeomDataHBSplines();


    ~GeomDataHBSplines();

    double  ComputeCoord(int dir, double param)
    { return  (origin[dir] + gridLEN[dir]*param);     }

    void  ComputeCoord(const myPoint& param, myPoint& geom)
    { 
      for(int ii=0; ii<DIM; ii++)
        geom[ii] = origin[ii] + gridLEN[ii]*param[ii];
      return;
    }

    double  ComputeParam(int dir, double val)
    { return  ((val-origin[dir])/gridLEN[dir]);     }

    void  ComputeParam(const myPoint& geom, myPoint& param)
    {
      for(int ii=0; ii<DIM; ii++)
        param[ii] = (geom[ii]-origin[ii])/gridLEN[ii];

      return;
    }

    void  SetDimension(int dim_)
    {  DIM = dim_;    }

    void  SetNDOF(int d)
    {  ndof = d;     }

    void  SetDegree(int ind, int deg)
    {   degree[ind] = deg;  }

    void  SetNGP(int ind, int nn)
    {    nGP[ind] = nn;     }

    void  SetNelem(int ind, int nn)
    {    nelem[ind] = nn;     }

    void  SetJacobian(int ind, double J1)
    {      Jacobian[ind] = J1;     }

    void  SetGridLength(int ind, double L1)
    {       gridLEN[ind] = L1;     }

    void  SetOrigin(int ind, double L1)
    {       origin[ind] = L1;     }

    int  GetNDOF()
    {       return  ndof;     }

    int GetDegree(int ind)
    {       return degree[ind];     }

    int GetNGP(int ind)
    {       return  nGP[ind];     }

    int GetNGBF(int ind)
    {       return  ngbf[ind];     }

    double  GetJacobian(int ind)
    {       return  Jacobian[ind];     }

    double  GetGridLength(int ind)
    {       return  gridLEN[ind];     }

    double  GetGridDX(int ind)
    {       return  dX[ind];     }

    double  GetJacobianFull()
    {       return  Jfull;     }

    int GetTotalNGP()
    {       return  totalNGP;     }

     void SetTimeIncrementType(int ttt)
     {  tis = ttt; }

     void SetRho(double ttt)
     {  rho = ttt; }

    void  printSelf();
    void  build();

    void  getBoundaryNormal1D(int side, myPoint& normal);
    void  getBoundaryNormal2D(int side, myPoint& normal);
    void  getBoundaryNormal3D(int side, myPoint& normal);

    void  setBoundaryGPs1D(int side, vector<double>& boundaryGPs1, vector<double>& boundaryGWs1);
    void  setBoundaryGPs2D(int side, vector<double>& boundaryGPs1, vector<double>& boundaryGWs1, vector<double>& boundaryGPs2, vector<double>& boundaryGWs2);
    void  setBoundaryGPs3D(int side, vector<double>& boundaryGPs1, vector<double>& boundaryGWs1, vector<double>& boundaryGPs2, vector<double>& boundaryGWs2, vector<double>& boundaryGPs3, vector<double>& boundaryGWs3);

    void  setBoundaryGPs1D(int side, GaussQuadrature& quadTemp);
    void  setBoundaryGPs2D(int side, GaussQuadrature& quadTemp);
    void  setBoundaryGPs3D(int side, GaussQuadrature& quadTemp);

    void  initialise(int size1=0, int size2=0, int size3=0, int size4=0);

    void  computeBasisFunctions1D(const myPoint& start, const myPoint& incr, const myPoint& param, VectorXd& N);

    void  computeBasisFunctions1D(const myPoint& start, const myPoint& incr, const myPoint& param, VectorXd& N, VectorXd& dN_dx);

    void  computeBasisFunctions1D(const myPoint& start, const myPoint& incr, const myPoint& param, VectorXd& N, VectorXd& dN_dx, VectorXd& d2N_dx2);

    void  computeBasisFunctions1D(const myPoint& start, const myPoint& incr, const myPoint& param, VectorXd& N, VectorXd& dN_dx, VectorXd& d2N_dx2, VectorXd& d3N_dx3);

    void  computeBasisFunctions2D(const myPoint& start, const myPoint& incr, const myPoint& param, VectorXd& N);

    void  computeBasisFunctions2D(const myPoint& start, const myPoint& incr, const myPoint& param, 
                              VectorXd& N, VectorXd& dN_dx, VectorXd& dN_dy);

    void  computeBasisFunctions2D(const myPoint& start, const myPoint& incr, const myPoint& param, 
                              VectorXd& N, VectorXd& dN_dx, VectorXd& dN_dy, VectorXd& d2N_dx2, VectorXd& d2N_dy2);

    void  computeBasisFunctions2D(const myPoint& start, const myPoint& incr, const myPoint& param, 
                              VectorXd& d3N_dx3, VectorXd& d3N_dy3, VectorXd& d3N_dxdy2, VectorXd& d3N_dx2dy);

    void  computeBasisFunctionsGhostPenalty2D(const myPoint& start, const myPoint& incr, const myPoint& param, 
                              VectorXd& N, VectorXd& dN_dx, VectorXd& dN_dy);

    void  computeBasisFunctions3D(const myPoint& start, const myPoint& incr, const myPoint& param, VectorXd& N);

    void  computeBasisFunctions3D(const myPoint& start, const myPoint& incr, const myPoint& param, VectorXd& N, VectorXd& dN_dx, VectorXd& dN_dy, VectorXd& dN_dz);

    void  computeBasisFunctions3D(const myPoint& start, const myPoint& incr, const myPoint& param, VectorXd& N, VectorXd& dN_dx, VectorXd& dN_dy, VectorXd& dN_dz,
                VectorXd& d2N_dx2, VectorXd& d2N_dy2, VectorXd& d2N_dz2);

    void  computeBasisFunctionsGhostPenalty3D(const myPoint& start, const myPoint& incr, const myPoint& param, 
                              VectorXd& N, VectorXd& dN_dx, VectorXd& dN_dy, VectorXd& dN_dz);

    void  reset();

    int  doIntersect2D(AABB& bbTemp, bool flag, vector<int>&  vecTemp,  vector<myPoint>& ptOut, vector<int>& domVec);
    
    int  doIntersect2Dfor3D(int sideTemp, double coord3, AABB& bbTemp, bool flag, vector<int>&  vecTemp,  vector<myPoint>& ptOut, vector<int>& domVec);

    int  doIntersect3D(AABB& bbTemp, bool flag, vector<int>&  vecTemp,  vector<myPoint>& ptOut, vector<int>& domVec);

    int  within(myPoint& pt);

};








#endif
