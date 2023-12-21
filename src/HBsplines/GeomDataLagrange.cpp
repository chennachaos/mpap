
#include "GeomDataLagrange.h"
#include "MpapTime.h"
#include "TimeFunction.h"
#include "BasisFunctionsLagrange.h"
#include "QuadratureUtil.h"


extern MpapTime mpapTime;
extern List<TimeFunction> timeFunction;



void GeomDataLagrange::build()
{
  //nGP = 1;
  
  //cout << " nGP " <<  nGP << endl;
  getGaussPoints1D(nGP, gausspoints1, gaussweights1);
  getGaussPoints1D(nGP, gausspoints2, gaussweights2);
  getGaussPoints1D(nGP, gausspoints3, gaussweights3);

  //printVector(gausspoints1);
  //printVector(gaussweights1);

  return;
}


void GeomDataLagrange::setNodalPositions(vector<myPoint>&  datatemp)
{
    int  ii=0, jj=0, ind=0;
    double  val=0.0;
  
    nNode = datatemp.size();

    NodePosOrig.resize(nNode);
    NodePosCur.resize(nNode);
    NodePosNew.resize(nNode);

    specValCur.resize(nNode);
    specValNew.resize(nNode);

    for(ii=0;ii<nNode;ii++)
    {
      for(jj=0;jj<DIM;jj++)
      {
        //cout << ii << '\t' << jj << '\t' << datatemp[ii][jj] << endl;
        val = datatemp[ii][jj];

        NodePosOrig[ii][jj] = val;
        NodePosCur[ii][jj]  = val;
        NodePosNew[ii][jj]  = val;

        specValCur[ii][jj]  = val;
        specValNew[ii][jj]  = val;
      }
    }

    if(DIM == 2)
    {
      for(ii=0;ii<nNode;ii++)
      {
        NodePosOrig[ii][jj] = 0.0;
        NodePosCur[ii][jj]  = 0.0;
        NodePosNew[ii][jj]  = 0.0;

        specValCur[ii][jj]  = 0.0;
        specValNew[ii][jj]  = 0.0;
      }
    }

    return;
}



void GeomDataLagrange::updateNodePositions(double* data)
{
    int  ii=0, jj=0, ind=0;
  
    nNode = NodePosCur.size();

    for(ii=0;ii<nNode;ii++)
    {
      ind = DIM*ii;
      for(jj=0;jj<DIM;jj++)
      {
        NodePosNew[ii][jj] = NodePosOrig[ii][jj] + data[ind+jj];
        NodePosCur[ii][jj] = NodePosOrig[ii][jj] + data[ind+jj];
      }
    }

    return;
}

void GeomDataLagrange::printSelf()
{
  //   cout << " Degree and Jacobian " << degree[0] << '\t' << Jfull << endl;
  return;
}



void GeomDataLagrange::reset()
{
  return;
}



void GeomDataLagrange::computeBasisFunctions1D(double uu, double *N, double *dN_dx)
{
    int  ROWS = 2, ii=0;
/*
    double** ders = new double*[ROWS], du_dx;

    for(ii=0;ii<ROWS;ii++)
       ders[ii] = new double[degree[0]+1];

    //HB_DersBasisFuns(degree[0], start, incr, uu, 1, ders);

    //  dx_du = J;
    //  du_dx = 1.0/dx_du

    //double  du_dx = 1.0/Jacobian[0];
    
    for(ii=0; ii<=degree[0]; ii++)
    {
        N[ii]     =  ders[0][ii];
        dN_dx[ii] =  ders[1][ii] * du_dx;
    }

    for(ii=0;ii<ROWS;ii++)
      delete [] ders[ii];

    delete [] ders;
*/
    return;
}




void GeomDataLagrange::computeBasisFunctions2D(int deg, double* param, double *N)
{
    int  ii=0, jj=0, count = deg + 1;

    vector<double>  N1(count), N2(count);

    Lagrange_BasisFuns1D(deg, param[0], &N1[0]);
    Lagrange_BasisFuns1D(deg, param[1], &N2[0]);

    count = 0;
    for(jj=0; jj<=deg; jj++)
    {
      for(ii=0; ii<=deg; ii++)
        N[count++]   =  N2[jj] *  N1[ii];
    }

    return;
}




void  GeomDataLagrange::computeBasisFunctions2D(bool flag, int type, int degree, double* param, vector<int>& nodeNums, double *N, double *dN_dx, double *dN_dy, double& Jac)
{
    int  ii=0, jj=0, count=0, nlbf=0;

    if(type == 1) // triangular elements
    {
      count = 3;
      nlbf = count;
    }
    else
    {
      count = degree + 1;
      nlbf = count*count;
    }

    vector<double>  N1(count), N2(count), dN1(count), dN2(count);
    vector<double>  dN_du1(nlbf), dN_du2(nlbf);
    double  xx, yy, detinv, B[2][2], Binv[2][2] ;

    if(type == 1) // triangular elements
      LagrangeBasisFunsTria(degree, param[0], param[1], &N[0], &dN_du1[0], &dN_du2[0]);
    else  // quad elements
    {
      Lagrange_BasisFuns1D(degree, param[0], &N1[0], &dN1[0]);
      Lagrange_BasisFuns1D(degree, param[1], &N2[0], &dN2[0]);

      count = 0;
      for(jj=0; jj<=degree; jj++)
      {
        for(ii=0; ii<=degree; ii++)
        {
          N[count]        =  N2[jj] *  N1[ii];
          dN_du1[count] =  N2[jj] * dN1[ii];
          dN_du2[count] = dN2[jj] *  N1[ii];
          count++;
        }
      }
    }

    B[0][0] = B[1][0] = B[0][1] = B[1][1] = 0.0;

    // Gradient of mapping from parameter space to physical space
    
    if(flag)
    {
       for(ii=0; ii<nlbf; ii++)
       {
          xx = NodePosCur[node_map_get_old[nodeNums[ii]]][0];
          yy = NodePosCur[node_map_get_old[nodeNums[ii]]][1];
          //cout << " coords ... " << xx << '\t' << yy << '\t' << dN_du1[ii] << '\t' << dN_du2[ii] << endl;

          B[0][0] +=  (xx * dN_du1[ii]) ;
          B[1][0] +=  (xx * dN_du2[ii]) ;
          B[0][1] +=  (yy * dN_du1[ii]) ;
          B[1][1] +=  (yy * dN_du2[ii]) ;
       }
    }
    else
    {
       for(ii=0; ii<nlbf; ii++)
       {
          xx = NodePosOrig[node_map_get_old[nodeNums[ii]]][0];
          yy = NodePosOrig[node_map_get_old[nodeNums[ii]]][1];

          //cout << " coords ... " << xx << '\t' << yy << '\t' << dN_du1[ii] << '\t' << dN_du2[ii] << endl;

          B[0][0] +=  (xx * dN_du1[ii]) ;
          B[1][0] +=  (xx * dN_du2[ii]) ;
          B[0][1] +=  (yy * dN_du1[ii]) ;
          B[1][1] +=  (yy * dN_du2[ii]) ;
       }
    }

    Jac  = B[0][0]*B[1][1] - B[0][1]*B[1][0];

    //printf("Bmat  \t%20.18f\t%20.18f\t%20.18f\t%20.18f\t%20.18f \n\n", B[0][0], B[0][1], B[1][0], B[1][1], Jac);

    detinv = 1.0/Jac ;

    Binv[0][0] =  B[1][1] * detinv;
    Binv[0][1] = -B[0][1] * detinv;
    Binv[1][0] = -B[1][0] * detinv;
    Binv[1][1] =  B[0][0] * detinv;

    // Compute derivatives of basis functions w.r.t physical coordinates
    for(ii=0; ii<nlbf; ii++)
    {
      dN_dx[ii] = dN_du1[ii] * Binv[0][0] + dN_du2[ii] * Binv[0][1];
      dN_dy[ii] = dN_du1[ii] * Binv[1][0] + dN_du2[ii] * Binv[1][1];
    }
  
    return;
}



void  GeomDataLagrange::computeDeformationGradient2D(bool flag, vector<int>& nodeNums, double* dN_dx, double* dN_dy, double* F, double& detF)
{
    /////////////////////////////////////////////
    //                                         //
    //    F(1,1) = F[0]   //   F(1,2) = F[1]   //
    //                                         //
    //    F(2,1) = F[2]   //   F(2,2) = F[3]   //
    //                                         //
    /////////////////////////////////////////////

    //   confflag = 1 --> coordiantes from current configuration
    //                    shape function derivatives w.r.t reference configuration

    //   confflag = 2 ==> coordiantes from reference configuration
    //                    shape function derivatives w.r.t current configuration

    double  xx=0.0, yy=0.0;

    F[0] = F[1] = F[2] = F[3] = 0.0;

    for(int ii=0; ii<nodeNums.size(); ii++)
    {
      xx = NodePosCur[node_map_get_old[nodeNums[ii]]][0];
      yy = NodePosCur[node_map_get_old[nodeNums[ii]]][1];

      //xx = NodePosCur[nodeNums[ii]][0];
      //yy = NodePosCur[nodeNums[ii]][1];
      //cout << xx << '\t' << yy << endl;

      F[0] += xx * dN_dx[ii];
      F[2] += xx * dN_dy[ii];
      F[1] += yy * dN_dx[ii];
      F[3] += yy * dN_dy[ii];
    }

    detF = F[0]*F[3] - F[1]*F[2];

    return;
}





void GeomDataLagrange::computeBasisFunctions3D(int deg, double* param, double *N)
{
    int  ii=0, jj=0, kk=0, count = deg + 1;

    vector<double>  N1(count), N2(count), N3(count);

    Lagrange_BasisFuns1D(deg, param[0], &N1[0]);
    Lagrange_BasisFuns1D(deg, param[1], &N2[0]);
    Lagrange_BasisFuns1D(deg, param[2], &N3[0]);

    count = 0;
    for(kk=0; kk<=deg; kk++)
    {
      for(jj=0; jj<=deg; jj++)
      {
        for(ii=0; ii<=deg; ii++)
          N[count++]   =  N3[kk] * N2[jj] * N1[ii];
      }
    }

    return;
}





void  GeomDataLagrange::computeBasisFunctions3D(bool flag, int type, int degree, double* param, vector<int>& nodeNums, double *N, double *dN_dx, double *dN_dy, double *dN_dz, double& Jac)
{
    int  ii=0, jj=0, kk=0, count=0, nlbf=0;

    if(type == 1) // tetrahedral elements
    {
      count = 4;
      nlbf = count;
    }
    else
    {
      count = degree + 1;
      nlbf = count*count*count;
    }

    vector<double>  N1(count), N2(count), N3(count), dN1(count), dN2(count), dN3(count);
    vector<double>  dN_du1(nlbf), dN_du2(nlbf), dN_du3(nlbf);
    double  xx=0.0, yy=0.0, zz=0.0;
    MatrixXd  B(3,3), Binv(3,3);

    if(type == 1) // tetrahedral elements
      LagrangeBasisFunsTet(degree, param[0], param[1], param[2], &N[0], &dN_du1[0], &dN_du2[0], &dN_du3[0]);
    else  // hex elements
    {
      Lagrange_BasisFuns1D(degree, param[0], &N1[0], &dN1[0]);
      Lagrange_BasisFuns1D(degree, param[1], &N2[0], &dN2[0]);
      Lagrange_BasisFuns1D(degree, param[2], &N3[0], &dN3[0]);

      count = 0;
      for(kk=0; kk<=degree; kk++)
      {
        for(jj=0; jj<=degree; jj++)
        {
          for(ii=0; ii<=degree; ii++)
          {
            N[count]       =  N3[kk] *  N2[jj] *  N1[ii];
            dN_du1[count]  =  N3[kk] *  N2[jj] * dN1[ii];
            dN_du2[count]  =  N3[kk] * dN2[jj] *  N1[ii];
            dN_du3[count]  = dN3[kk] *  N2[jj] *  N1[ii];

            count++;
          }
        }
      }
    }

    // Gradient of mapping from parameter space to physical space
    B.setZero();

    if(flag)
    {
       for(ii=0; ii<nlbf; ii++)
       {
          xx = NodePosCur[node_map_get_old[nodeNums[ii]]][0];
          yy = NodePosCur[node_map_get_old[nodeNums[ii]]][1];
          zz = NodePosCur[node_map_get_old[nodeNums[ii]]][2];

          //xx = NodePosCur[nodeNums[ii]][0];
          //yy = NodePosCur[nodeNums[ii]][1];
          //zz = NodePosCur[nodeNums[ii]][2];
          //cout << " coords ... " << xx << '\t' << yy << '\t' << dN_du1[ii] << '\t' << dN_du2[ii] << endl;

          B(0,0) +=  (xx * dN_du1[ii]) ;
          B(1,0) +=  (xx * dN_du2[ii]) ;
          B(2,0) +=  (xx * dN_du3[ii]) ;

          B(0,1) +=  (yy * dN_du1[ii]) ;
          B(1,1) +=  (yy * dN_du2[ii]) ;
          B(2,1) +=  (yy * dN_du3[ii]) ;

          B(0,2) +=  (zz * dN_du1[ii]) ;
          B(1,2) +=  (zz * dN_du2[ii]) ;
          B(2,2) +=  (zz * dN_du3[ii]) ;
      }
    }
    else
    {
      for(ii=0; ii<nlbf; ii++)
      {
          xx = NodePosOrig[node_map_get_old[nodeNums[ii]]][0];
          yy = NodePosOrig[node_map_get_old[nodeNums[ii]]][1];
          zz = NodePosOrig[node_map_get_old[nodeNums[ii]]][2];

          //cout << " coords ... " << xx << '\t' << yy << '\t' << dN_du1[ii] << '\t' << dN_du2[ii] << endl;

          B(0,0) +=  (xx * dN_du1[ii]) ;
          B(1,0) +=  (xx * dN_du2[ii]) ;
          B(2,0) +=  (xx * dN_du3[ii]) ;

          B(0,1) +=  (yy * dN_du1[ii]) ;
          B(1,1) +=  (yy * dN_du2[ii]) ;
          B(2,1) +=  (yy * dN_du3[ii]) ;

          B(0,2) +=  (zz * dN_du1[ii]) ;
          B(1,2) +=  (zz * dN_du2[ii]) ;
          B(2,2) +=  (zz * dN_du3[ii]) ;
      }
    }

    Jac  = B.determinant();
    Binv = B.inverse();

    // Compute derivatives of basis functions w.r.t physical coordinates
    for(ii=0; ii<nlbf; ii++)
    {
      dN_dx[ii] = dN_du1[ii] * Binv(0,0) + dN_du2[ii] * Binv(0,1) + dN_du3[ii] * Binv(0,2);
      dN_dy[ii] = dN_du1[ii] * Binv(1,0) + dN_du2[ii] * Binv(1,1) + dN_du3[ii] * Binv(1,2);
      dN_dz[ii] = dN_du1[ii] * Binv(2,0) + dN_du2[ii] * Binv(2,1) + dN_du3[ii] * Binv(2,2);
    }

    return;
}




void  GeomDataLagrange::computeDeformationGradient3D(bool flag, vector<int>& nodeNums, double* dN_dx, double* dN_dy, double* dN_dz, double* F, double& detF)
{
    double  xx=0.0, yy=0.0, zz=0.0;

    F[0] = F[1] = F[2] = 0.0;
    F[3] = F[4] = F[5] = 0.0;
    F[6] = F[7] = F[8] = 0.0;

    for(int ii=0; ii<nodeNums.size(); ii++)
    {
        xx = NodePosCur[node_map_get_old[nodeNums[ii]]][0];
        yy = NodePosCur[node_map_get_old[nodeNums[ii]]][1];
        zz = NodePosCur[node_map_get_old[nodeNums[ii]]][2];

        //cout << xx << '\t' << yy << endl;

        F[0] += xx * dN_dx[ii];
        F[3] += xx * dN_dy[ii];
        F[6] += xx * dN_dz[ii];

        F[1] += yy * dN_dx[ii];
        F[4] += yy * dN_dy[ii];
        F[7] += yy * dN_dz[ii];

        F[2] += zz * dN_dx[ii];
        F[5] += zz * dN_dy[ii];
        F[8] += zz * dN_dz[ii];
    }

    detF = F[0]*(F[4]*F[8] - F[5]*F[7]) - F[3]*(F[1]*F[8] - F[2]*F[7]) + F[6]*(F[1]*F[5] - F[2]*F[4]);

    return;
}





