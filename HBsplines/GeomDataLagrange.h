
#ifndef incl_GeomDataLagrange_h
#define incl_GeomDataLagrange_h

#include "util.h"
#include <vector>

using std::vector;


class  GeomDataLagrange
{
  private:
    int  DIM, nGP, nNode, ndof;

  public:

    vector<double>  gausspoints1, gaussweights1, gausspoints2, gaussweights2, gausspoints3, gaussweights3;

    vector<myPoint>  gaussps, gaussws;

    vector<myPoint>  NodePosOrig, NodePosNew, NodePosCur, NodePosPrev, NodePosPrevCur;
    vector<myPoint>  specValNew, specValCur, acceNew, acceCur;

    vector<int>  node_map_new_to_old, node_map_old_to_new;
    vector<int>  assy4r;

    GeomDataLagrange()  {}

    ~GeomDataLagrange() {}

    void  setDimension(int dd)
    {      DIM = dd;    }

    void  setNumNode(int dd)
    {      nNode = dd;    }

    void  setNdof(int dd)
    {      ndof = dd;    }

    void  setNGP(int nn)
    {      nGP = nn;    }

    int  getNumNode()
    {  return      nNode;    }

    int  getNdof()
    {  return  ndof;     }

    int getNGP(int ind)
    {  return  nGP;     }

    void  printSelf();

    void  build();

    void  reset();

    void  setNodalPositions(vector<myPoint>&  datatemp);

    void  updateNodePositions(double* );

    void  computeBasisFunctions1D(double uu, double *N, double *dN_dx);

    void  computeBasisFunctions2D(int deg, double* param, double *N);

    void  computeBasisFunctions2D(bool flag, int type, int deg, double* param, vector<int>& nodeNum, double *N, double *dN_dx, double *dN_dy, double& Jac);

    void  computeDeformationGradient2D(bool flag, vector<int>& nodeNum, double* dN_dx, double* dN_dy, double* F, double& detF);

    void  computeBasisFunctions3D(int deg, double* param, double *N);

    void  computeBasisFunctions3D(bool flag, int type, int deg, double* param, vector<int>& nodeNum, double *N, double *dN_dx, double *dN_dy, double* dN_dz, double& Jac);

    void  computeDeformationGradient3D(bool flag, vector<int>& nodeNum, double* dN_dx, double* dN_dy, double* dN_dz, double* F, double& detF);

};








#endif
