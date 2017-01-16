#include <Eigen/Dense>
#include "NurbsSURFACE.h"
#include "DataBlockTemplate.h"
#include "Plot.h"
#include "PlotVTK.h"
#include <iomanip>

#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataMapper.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkCellArray.h>
#include <vtkXMLPolyDataWriter.h>  
#include <vtkPolygon.h>  
#include <vtkActor2D.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataSetMapper.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkProperty.h>
#include <vtkFloatArray.h>
#include <vtkLine.h>
#include <vtkQuad.h>
#include <vtkScalarBarActor.h>
#include <vtkTextProperty.h>
#include <vtkLookupTable.h>



using namespace std;
using namespace Eigen;

extern Plot plot;
extern PlotVTK plotvtk;
extern MpapTime mpapTime;








double  NurbsSURFACE::computeValue(int dof, double u, double v)
{
    vector<double>  NN(nlbf);

    return computeValueAndShanpeFns(dof, u, v, &NN[0]);
}





double  NurbsSURFACE::computeValueAndShanpeFns(int dof, double u, double v, double* NN)
{
  double  val=0.0, temp;

    double  Nu[nlbf1], Nv[nlbf2];

    int ni, nj, ind, ii, jj, count, loc_num;

    ni = FindSpan(&(U[0]), U.n, p, u) - p ;
    nj = FindSpan(&(V[0]), V.n, q, v) - q;

    BasisFuns(&(U[0]), U.n, p, u, Nu) ;
    BasisFuns(&(V[0]), V.n, q, v, Nv) ;

    double sum_tot = 0.0, wt;
    loc_num = 0;

    for(jj=0; jj<=q; jj++)
    {
      temp = nj+jj;

      for(ii=0; ii<=p; ii++)
      {
        wt = Pw[ni+ii][temp].w;

        NN[loc_num] = Nv[jj] * Nu[ii] * wt;

        sum_tot  +=  NN[loc_num];

        loc_num++;
      }
    }

    loc_num = 0;
    val = 0.0;
    dof -= 1;

    for(jj=0; jj<=q; jj++)
    {
      ind = ni + ngbf1 * (nj+jj);
      for(ii=0; ii<=p; ii++)
      {
        NN[loc_num] = NN[loc_num]/sum_tot;

        val += NN[loc_num] * Values[dof][ind+ii];

        loc_num++;
      }
    }

  return val;
}




void NurbsSURFACE::ShapeFunctions(double u, double v , double* NN)
{

    double Nu[nlbf1], Nv[nlbf2];

    BasisFuns(&(U[0]), U.n, p, u, Nu);
    BasisFuns(&(V[0]), V.n, q, v, Nv);

    int loc_num = 0, ii, jj, temp, ni, nj;

    ni = FindSpan(&(U[0]), U.n, p, u) - p;
    nj = FindSpan(&(V[0]), V.n, q, v) - q;

    double sum_tot = 0.0, wt;
    loc_num = 0;

    for(jj=0; jj<=q; jj++)
    {
      temp = nj+jj;

      for(ii=0; ii<=p; ii++)
      {
        wt = Pw[ni+ii][temp].w;

        NN[loc_num] = Nu[ii] * Nv[jj] * wt;

        sum_tot  +=  NN[loc_num];

        loc_num++;
      }
    }

    for(loc_num = 0;loc_num<nlbf; loc_num++)
      NN[loc_num] = NN[loc_num]/sum_tot;

  return;
}




void NurbsSURFACE::ShapeFunDerivatives5(int index, double* params, double* dN_dxi)
{
    int   ii, jj, loc_num, ROWS = 2;
   
    double** ders1 = new double*[ROWS];
    double** ders2 = new double*[ROWS];

    for(ii=0;ii<ROWS;ii++)
    {
       ders1[ii] = new double[nlbf1];
       ders2[ii] = new double[nlbf2];
    }

    DersBasisFuns(&(U[0]), U.n, p, params[0], 1, ders1 );
    DersBasisFuns(&(V[0]), V.n, q, params[1], 1, ders2 );

    if(index == 1)
    {
       loc_num = 0;
       for(jj=0; jj<=q; jj++)
       {
          for(ii=0; ii<=p; ii++)
          {
             dN_dxi[loc_num]  = ders1[1][ii] * ders2[0][jj];
             loc_num++;
          }
       }
    }
    else
    {
       loc_num = 0;
       for(jj=0; jj<=q; jj++)
       {
          for(ii=0; ii<=p; ii++)
          {
             dN_dxi[loc_num]  = ders1[0][ii] * ders2[1][jj];
             loc_num++;
          }
       }
    }


  for(ii=0;ii<ROWS;ii++)
  {
    delete [] ders1[ii];
    delete [] ders2[ii];
  }
  delete [] ders1;
  delete [] ders2;

  return;
}





void NurbsSURFACE::ShapeFunDerivatives(int* startindex, double* params, double* N, double* dN_dx, double* dN_dy, double& Jac)
{
    int   temp, ni, nj, ii, jj, kk, loc_num, ind1, ind2, ind3, ROWS = 2;

    double  det, detinv, dR_du[nlbf][2], B[2][2], Binv[2][2], sum_tot=0.0, wt, sum1, sum2 ;

    EPOINT *EP;

    ni = startindex[0];
    nj = startindex[1];
   
    double** ders1 = new double*[ROWS];
    double** ders2 = new double*[ROWS];

    for(ii=0;ii<ROWS;ii++)
    {
       ders1[ii] = new double[nlbf1];
       ders2[ii] = new double[nlbf2];
    }

    DersBasisFuns(&(U[0]), U.n, p, params[0], 1, ders1 );
    DersBasisFuns(&(V[0]), V.n, q, params[1], 1, ders2 );

    // Calculate Bivariate B-Spline basis functions and their derivatives w.r.t parametetric coordinates

    sum1 = sum2 = 0.0;
    sum_tot = 0.0;
    loc_num = 0;
    for(jj=0; jj<=q; jj++)
    {
        temp = nj+jj;
        for(ii=0; ii<=p; ii++)
        {
          wt = Pw[ni+ii][temp].w;

          N[loc_num] = ders1[0][ii] * ders2[0][jj] * wt;
          
          sum_tot  +=  N[loc_num];
       
          dR_du[loc_num][0] = ders1[1][ii] * ders2[0][jj] * wt;
          dR_du[loc_num][1] = ders1[0][ii] * ders2[1][jj] * wt;
          
          sum1 += dR_du[loc_num][0];
          sum2 += dR_du[loc_num][1];

          loc_num++;
        }
    }
    
    det = sum_tot*sum_tot;

    for(loc_num = 0;loc_num<nlbf; loc_num++)
    {
      dR_du[loc_num][0] = (dR_du[loc_num][0]*sum_tot - N[loc_num]*sum1)/det;
      dR_du[loc_num][1] = (dR_du[loc_num][1]*sum_tot - N[loc_num]*sum2)/det;

      N[loc_num] = N[loc_num]/sum_tot;
    }


    B[0][0] = B[1][0] = B[0][1] = B[1][1] = 0.0;

    // Gradient of mapping from parameter space to physical space
    loc_num = 0;

    for(jj=0; jj<=q; jj++)
    {
       temp = nj+jj;
       for(ii=0; ii<=p; ii++)
       {
          EP = &(PP[ni+ii][temp] );

          B[0][0] +=  (EP->x * dR_du[loc_num][0]) ;
          B[1][0] +=  (EP->x * dR_du[loc_num][1]) ;
          B[0][1] +=  (EP->y * dR_du[loc_num][0]) ;
          B[1][1] +=  (EP->y * dR_du[loc_num][1]) ;

          loc_num++;
       }
    }

  Jac = B[0][0]*B[1][1] - B[0][1]*B[1][0];
  
  detinv = 1.0/Jac ;

  Binv[0][0] =  B[1][1] * detinv;
  Binv[0][1] = -B[0][1] * detinv;
  Binv[1][0] = -B[1][0] * detinv;
  Binv[1][1] =  B[0][0] * detinv;

  // Compute derivatives of basis functions w.r.t physical coordinates
  for(loc_num = 0;loc_num<nlbf; loc_num++)
  {
    dN_dx[loc_num] = dR_du[loc_num][0] * Binv[0][0] + dR_du[loc_num][1] * Binv[0][1];
    dN_dy[loc_num] = dR_du[loc_num][0] * Binv[1][0] + dR_du[loc_num][1] * Binv[1][1];
  }
  
  EP = NULL;
  for(ii=0;ii<ROWS;ii++)
  {
    delete [] ders1[ii];
    delete [] ders2[ii];
  }
  delete [] ders1;
  delete [] ders2;

  return;
}




/*
void NurbsSURFACE::ShapeFunDerivatives2(int* startindex, double* params, double* N, double* dN_dx, double* dN_dy, double* d2N_dx2, double* d2N_dy2, double& Jac)
{
    int   temp, ni, nj, ii, jj, kk, loc_num, ind1, ind2, ind3, ROWS = 3;
   
    double  dNdu[nlbf], dMdv[nlbf], d2x[3], d2y[3], d2N[nlbf], d2M[nlbf], d1N1[nlbf];

    double  dJdx, dJdy, d2udx2, d2udy2, d2vdx2, d2vdy2, fact;
    
    MatrixXd  B(2,2), Binv(2,2);
   
    ni = startindex[0];
    nj = startindex[1];
   
    double** ders1 = new double*[ROWS];
    double** ders2 = new double*[ROWS];

    for(ii=0;ii<ROWS;ii++)
    {
       ders1[ii] = new double[nlbf1];
       ders2[ii] = new double[nlbf2];
    }

    DersBasisFuns(&(U[0]), U.n, p, params[0], 2, ders1 );
    DersBasisFuns(&(V[0]), V.n, q, params[1], 2, ders2 );

    // Calculate Bivariate B-Spline basis functions and their derivatives w.r.t parametetric coordinates

    loc_num = 0;
    for(jj=0; jj<=q; jj++)
    {
       for(ii=0; ii<=p; ii++)
       {
          N[loc_num]     =  ders1[0][ii] * ders2[0][jj];
          dNdu[loc_num]  =  ders1[1][ii] * ders2[0][jj];
          dMdv[loc_num]  =  ders1[0][ii] * ders2[1][jj];
          
          d2N[loc_num]   =  ders1[2][ii] * ders2[0][jj];
          d2M[loc_num]   =  ders1[0][ii] * ders2[2][jj];
          d1N1[loc_num]  =  ders1[1][ii] * ders2[1][jj];
          
          loc_num++;
       }
    }

    B.setZero();
    
    d2x[0] = d2x[1] = d2x[2] = 0.0;
    d2y[0] = d2y[1] = d2y[2] = 0.0;

    // Gradient of mapping from parameter space to physical space
    loc_num = 0;
    EPOINT *EP;

    for(jj=0; jj<=q; jj++)
    {
       temp = nj+jj;
       for(ii=0; ii<=p; ii++)
       {
          EP = &(PP[ni+ii][temp] );

          B(0,0) +=  (EP->x * dNdu[loc_num]) ;
          B(1,0) +=  (EP->x * dMdv[loc_num]) ;
          B(0,1) +=  (EP->y * dNdu[loc_num]) ;
          B(1,1) +=  (EP->y * dMdv[loc_num]) ;
          
          d2x[0]  +=  (EP->x * d2N[loc_num]) ;
          d2x[1]  +=  (EP->x * d1N1[loc_num]) ;
          d2x[2]  +=  (EP->x * d2M[loc_num]) ;

          d2y[0]  +=  (EP->y * d2N[loc_num]) ;
          d2y[1]  +=  (EP->y * d1N1[loc_num]) ;
          d2y[2]  +=  (EP->y * d2M[loc_num]) ;

          loc_num++;
       }
    }

    //printf(" \t %14.8f\t %14.8f\t%14.8f\t%14.8f\n", B(0,0),B(0,1),B(1,0),B(1,1));
    //printf(" \t %14.8f\t %14.8f\t%14.8f\t %14.8f\t%14.8f\t%14.8f\n\n", d2x[0], d2x[1], d2x[2], d2y[0], d2y[1], d2y[2]);

  Jac = B.determinant();

  fact = Jac*Jac;
  
  Binv = B.inverse();

  // Compute derivatives of basis functions w.r.t physical coordinates
  for(loc_num = 0;loc_num<nlbf; loc_num++)
  {
    dN_dx[loc_num] = dNdu[loc_num] * Binv(0,0) + dMdv[loc_num] * Binv(0,1);
    dN_dy[loc_num] = dNdu[loc_num] * Binv(1,0) + dMdv[loc_num] * Binv(1,1);
    
    dJdx  =  (d2x[0]*B(1,1) - d2x[1]*B(0,1))*B(1,1) + (d2y[1]*B(1,1) - d2y[2]*B(0,1))*B(0,0);
    dJdx -= ((d2x[1]*B(1,1) - d2x[2]*B(0,1))*B(0,1) + (d2y[0]*B(1,1) - d2y[1]*B(0,1))*B(1,0));
    
    dJdx /= Jac;

    dJdy  =  (-d2x[0]*B(1,0) + d2x[1]*B(0,0))*B(1,1) + (-d2y[1]*B(1,0) + d2y[2]*B(0,0))*B(0,0);
    dJdy += (-(d2x[1]*B(1,0) + d2x[2]*B(0,0))*B(0,1) - (-d2y[0]*B(1,0) + d2y[1]*B(0,0))*B(1,0));
    
    dJdy /= Jac;

    d2udx2 = ((d2y[1]- dJdx)*B(1,1) - d2y[2]*B(0,1))/fact ;

    d2vdx2 = (-d2y[0]*B(1,1) + (d2y[1] + dJdx)*B(0,1))/fact ;
    
    d2udy2 = ((d2x[1]+dJdy)*B(1,0) - d2x[2]*B(0,0))/fact ;

    d2vdy2 = (-d2x[0]*B(1,0) + (d2x[1]-dJdy)*B(0,0))/fact ;

    //printf(" \t %14.8f\t %14.8f\t%14.8f\t %14.8f\t%14.8f\t%14.8f\n", dJdx, dJdy, d2udx2, d2udy2, d2vdx2, d2vdy2);

    d2N_dx2[loc_num]  =  dNdu[loc_num] * d2udx2 + d2N[loc_num] * Binv(0,0) * Binv(0,0) + 2.0*d1N1[loc_num] * Binv(0,0) * Binv(0,1) ;
    d2N_dx2[loc_num] += (dMdv[loc_num] * d2vdx2 + d2M[loc_num] * Binv(0,1) * Binv(0,1));
    
    d2N_dy2[loc_num]  =  dNdu[loc_num] * d2udy2 + d2N[loc_num] * Binv(1,0) * Binv(1,0) + 2.0*d1N1[loc_num] * Binv(1,1) * Binv(1,0) ;
    d2N_dy2[loc_num] += (dMdv[loc_num] * d2vdy2 + d2M[loc_num] * Binv(1,1) * Binv(1,1));
  }
  
//  cout << " UUUUUUUUUUUUU " << endl;

  EP = NULL;

  for(ii=0;ii<ROWS;ii++)
  {
    delete [] ders1[ii];
    delete [] ders2[ii];
  }
  delete [] ders1;
  delete [] ders2;

  return;
}
*/



void NurbsSURFACE::ShapeFunDerivatives2(int* startindex, double* params, double* N, double* dN_dx, double* dN_dy, double* d2N_dx2, double* d2N_dy2, double* d2N_dxy, double* d2N_dyx, double& Jac)
{

    int   temp, ni, nj, ii, jj, kk, loc_num, ind1, ind2, ind3, ROWS = 3;
   
    double  dNdu[nlbf], dNdv[nlbf], d2x[3], d2y[3], d2Ndu2[nlbf], d2Ndv2[nlbf], d2Ndudv[nlbf];

    MatrixXd  B(2,2), Binv(2,2), BinvT(2,2), Hx(2,2), Hy(2,2), HNx(2,2), HNxi(2,2), C(2,2);

    ni = startindex[0];
    nj = startindex[1];

    double** ders1 = new double*[ROWS];
    double** ders2 = new double*[ROWS];

    for(ii=0;ii<ROWS;ii++)
    {
       ders1[ii] = new double[nlbf1];
       ders2[ii] = new double[nlbf2];
    }

    DersBasisFuns(&(U[0]), U.n, p, params[0], 2, ders1 );
    DersBasisFuns(&(V[0]), V.n, q, params[1], 2, ders2 );

    // Calculate Bivariate B-Spline basis functions and their derivatives w.r.t parametetric coordinates
    //cout << " ppppppppppppppp " << endl;


    double  w=0.0, dwdu=0.0, dwdv=0.0, d2wdu2=0.0, d2wdv2=0.0, d2wdudv=0.0, wt;

    loc_num = 0;
    for(jj=0; jj<=q; jj++)
    {
        temp = nj+jj;
        for(ii=0; ii<=p; ii++)
        {
          wt = Pw[ni+ii][temp].w;

          N[loc_num]     =  ders2[0][jj] * ders1[0][ii] * wt;
          dNdu[loc_num]  =  ders2[0][jj] * ders1[1][ii] * wt;
          dNdv[loc_num]  =  ders2[1][jj] * ders1[0][ii] * wt;
          
          d2Ndu2[loc_num]   =  ders2[0][jj] * ders1[2][ii] * wt;
          d2Ndv2[loc_num]   =  ders2[2][jj] * ders1[0][ii] * wt;
          d2Ndudv[loc_num]  =  ders2[1][jj] * ders1[1][ii] * wt;

          w       += N[loc_num];
          dwdu    += dNdu[loc_num];
          dwdv    += dNdv[loc_num];
          d2wdu2  += d2Ndu2[loc_num];
          d2wdv2  += d2Ndv2[loc_num];
          d2wdudv += d2Ndudv[loc_num];

          loc_num++;
        }
    }

    double  w2 = w*w;
    double  w3 = w2*w;

    for(ii=0; ii<nlbf; ii++)
    {
      dNdu[ii] = (dNdu[ii]*w - N[ii]*dwdu)/w2;
      dNdv[ii] = (dNdv[ii]*w - N[ii]*dwdv)/w2;
      
      d2Ndu2[ii] = -2.0*(w*dNdu[ii] - dwdu*N[ii])*dwdu/w3 + (w*d2Ndu2[ii] - d2wdu2*N[ii])/w2;

      d2Ndv2[ii] = -2.0*(w*dNdv[ii] - dwdv*N[ii])*dwdv/w3 + (w*d2Ndv2[ii] - d2wdv2*N[ii])/w2;

      d2Ndudv[ii] = d2Ndudv[ii]/w - d2wdudv*N[ii]/w2 - dNdv[ii]*dwdu/w2 + 2.0*dwdu*dwdv*N[ii]/w3 - dwdv*dNdu[ii]/w2 ;

      N[ii] = N[ii]/w;
    }


    //cout << " ppppppppppppppp " << endl;

    B.setZero();    Hx.setZero();    Hy.setZero();

    // Gradient of mapping from parameter space to physical space
    loc_num = 0;
    EPOINT *EP;

    for(jj=0; jj<=q; jj++)
    {
       temp = nj+jj;
       for(ii=0; ii<=p; ii++)
       {
          EP = &(PP[ni+ii][temp] );

          B(0,0) +=  (EP->x * dNdu[loc_num]) ;
          B(1,0) +=  (EP->x * dNdv[loc_num]) ;
          B(0,1) +=  (EP->y * dNdu[loc_num]) ;
          B(1,1) +=  (EP->y * dNdv[loc_num]) ;
          
          Hx(0,0) +=  (EP->x * d2Ndu2[loc_num]) ;
          Hx(0,1) +=  (EP->x * d2Ndudv[loc_num]) ;
          Hx(1,1) +=  (EP->x * d2Ndv2[loc_num]) ;

          Hy(0,0)  +=  (EP->y * d2Ndu2[loc_num]) ;
          Hy(0,1)  +=  (EP->y * d2Ndudv[loc_num]) ;
          Hy(1,1)  +=  (EP->y * d2Ndv2[loc_num]) ;

          loc_num++;
       }
    }
    Hx(1,0) = Hx(0,1);
    Hy(1,0) = Hy(0,1);

    //cout << " qqqqqqqqqqqqqq " << endl;

  Jac = B.determinant();

  Binv = B.inverse();
  BinvT = Binv.transpose();

  /*
  cout << " Jac...." << B(0,0) << '\t' << B(0,1) << '\t' << B(1,0) << '\t' << B(1,1) << '\t' << Jac << endl;
  cout << endl;cout << endl;
  cout << " Hx....." << Hx(0,0) << '\t' << Hx(0,1) << '\t' << Hx(1,0) << '\t' << Hx(1,1) << endl;
  cout << endl;cout << endl;
  cout << " Hy....." << Hy(0,0) << '\t' << Hy(0,1) << '\t' << Hy(1,0) << '\t' << Hy(1,1) << endl;
  cout << endl;cout << endl;
  cout << " Binv....." << Binv(0,0) << '\t' << Binv(0,1) << '\t' << Binv(1,0) << '\t' << Binv(1,1) << endl;
  cout << endl;cout << endl;
  cout << " BinvT....." << BinvT(0,0) << '\t' << BinvT(0,1) << '\t' << BinvT(1,0) << '\t' << BinvT(1,1) << endl;
  cout << endl;cout << endl;
  cout << endl;cout << endl;
  */

  // Compute derivatives of basis functions w.r.t physical coordinates
  for(ii = 0;ii<nlbf;ii++)
  {
    dN_dx[ii] = Binv(0,0) * dNdu[ii] + Binv(0,1) * dNdv[ii] ;
    dN_dy[ii] = Binv(1,0) * dNdu[ii] + Binv(1,1) * dNdv[ii] ;

    HNxi(0,0) = d2Ndu2[ii];
    HNxi(0,1) = d2Ndudv[ii];
    HNxi(1,0) = d2Ndudv[ii];
    HNxi(1,1) = d2Ndv2[ii];

    C = Hx*dN_dx[ii] + Hy*dN_dy[ii];

    HNx = Binv * (HNxi - C) * BinvT;

    d2N_dx2[ii]  =  HNx(0,0);
    d2N_dxy[ii]  =  HNx(0,1);
    d2N_dyx[ii]  =  HNx(1,0);
    d2N_dy2[ii]  =  HNx(1,1);

    //printf("HNx..... \t%12.6f \t %12.6f \t %12.6f \t %12.6f \n", HNx(0,0), HNx(0,1), HNx(1,0), HNx(1,1));
  }
  //cout << endl;cout << endl;

    //cout << " qqqqqqqqqqqqqq " << endl;

  EP = NULL;

  for(ii=0;ii<ROWS;ii++)
  {
    delete [] ders1[ii];
    delete [] ders2[ii];
  }
  delete [] ders1;
  delete [] ders2;

  return;
}



/*
void NurbsSURFACE::ShapeFunDerivatives4(int* startindex, double* params, double* d2N_dx2, double* d2N_dy2, double* d2N_dxy, double* d2N_dyx, double& Jac)
{
    int   temp, ni, nj, ii, jj, kk, loc_num, ind1, ind2, ind3, ROWS = 3;
   
    double  dNdu[nlbf], dMdv[nlbf], d2x[3], d2y[3], d2N[nlbf], d2M[nlbf], d1N1[nlbf];
    double  N[nlbf], dN_dx[nlbf], dN_dy[nlbf];
    
    MatrixXd  B(2,2), Binv(2,2), BinvT(2,2), Hx(2,2), Hy(2,2), HNx(2,2), HNxi(2,2), C(2,2);
   
    ni = startindex[0];
    nj = startindex[1];
   
    double** ders1 = new double*[ROWS];
    double** ders2 = new double*[ROWS];

    for(ii=0;ii<ROWS;ii++)
    {
       ders1[ii] = new double[nlbf1];
       ders2[ii] = new double[nlbf2];
    }

    DersBasisFuns(&(U[0]), U.n, p, params[0], 2, ders1 );
    DersBasisFuns(&(V[0]), V.n, q, params[1], 2, ders2 );

    // Calculate Bivariate B-Spline basis functions and their derivatives w.r.t parametetric coordinates

    loc_num = 0;
    for(jj=0; jj<=q; jj++)
    {
       for(ii=0; ii<=p; ii++)
       {
          N[loc_num]     =  ders1[0][ii] * ders2[0][jj];
          dNdu[loc_num]  =  ders1[1][ii] * ders2[0][jj];
          dMdv[loc_num]  =  ders1[0][ii] * ders2[1][jj];
          
          d2N[loc_num]   =  ders1[2][ii] * ders2[0][jj];
          d2M[loc_num]   =  ders1[0][ii] * ders2[2][jj];
          d1N1[loc_num]  =  ders1[1][ii] * ders2[1][jj];
          
          loc_num++;
       }
    }

    B.setZero();    Hx.setZero();    Hy.setZero();

    // Gradient of mapping from parameter space to physical space
    loc_num = 0;
    EPOINT *EP;

    for(jj=0; jj<=q; jj++)
    {
       temp = nj+jj;
       for(ii=0; ii<=p; ii++)
       {
          EP = &(PP[ni+ii][temp] );

          B(0,0) +=  (EP->x * dNdu[loc_num]) ;
          B(1,0) +=  (EP->x * dMdv[loc_num]) ;
          B(0,1) +=  (EP->y * dNdu[loc_num]) ;
          B(1,1) +=  (EP->y * dMdv[loc_num]) ;
          
          Hx(0,0) +=  (EP->x * d2N[loc_num]) ;
          Hx(0,1) +=  (EP->x * d1N1[loc_num]) ;
          Hx(1,1)  +=  (EP->x * d2M[loc_num]) ;

          Hy(0,0)  +=  (EP->y * d2N[loc_num]) ;
          Hy(0,1)  +=  (EP->y * d1N1[loc_num]) ;
          Hy(1,1)  +=  (EP->y * d2M[loc_num]) ;

          loc_num++;
       }
    }
    Hx(1,0) = Hx(0,1);
    Hy(1,0) = Hy(0,1);

  Jac = B.determinant();
  
  Binv = B.inverse();
  BinvT = Binv.transpose();

  // Compute derivatives of basis functions w.r.t physical coordinates
  for(ii=0;ii<nlbf;ii++)
  {
    dN_dx[ii] = dNdu[ii] * Binv(0,0) + dMdv[ii] * Binv(0,1);
    dN_dy[ii] = dNdu[ii] * Binv(1,0) + dMdv[ii] * Binv(1,1);

    HNxi(0,0) = d2N[ii];
    HNxi(0,1) = d1N1[ii];
    HNxi(1,0) = d1N1[ii];
    HNxi(1,1) = d2M[ii];

    C = Hx*dN_dx[ii] + Hy*dN_dy[ii];

    HNx = Binv * (HNxi - C) * BinvT;

    d2N_dx2[ii]  =  HNx(0,0);
    d2N_dxy[ii]  =  HNx(0,1);
    d2N_dyx[ii]  =  HNx(1,0);
    d2N_dy2[ii]  =  HNx(1,1);
  }
  
  EP = NULL;

  for(ii=0;ii<ROWS;ii++)
  {
    delete [] ders1[ii];
    delete [] ders2[ii];
  }
  delete [] ders1;
  delete [] ders2;

  return;
}
*/




double  NurbsSURFACE::computeValue2(int ni, int nj, double* NN)
{
      int  ii, jj, loc_num, temp;
      
      double  val=0.0;

      EPOINT *EP;

         loc_num = 0;
         for(jj=0; jj<=q; jj++)
         {
            temp = nj+jj;
            for(ii=0; ii<=p; ii++)
            {
               EP = &(PP[ni+ii][temp] );

               val += EP->z * NN[loc_num++];
            }
         }

  return val;
}




//
void NurbsSURFACE::deformationGradient(int ni, int nj, bool confflag, double* dN_dx, double* dN_dy, double* F, double& detF)
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

      int  ii, jj, loc_num, temp;

      EPOINT *EP;

      F[0] = F[1] = F[2] = F[3] = 0.0;

      if(confflag) // calculation of F
      {
         loc_num = 0;
         for(jj=0; jj<=q; jj++)
         {
            temp = nj+jj;
            for(ii=0; ii<=p; ii++)
            {
               EP = &(PP[ni+ii][temp] );

               F[0] += EP->x * dN_dx[loc_num];
               F[2] += EP->x * dN_dy[loc_num];
               F[1] += EP->y * dN_dx[loc_num];
               F[3] += EP->y * dN_dy[loc_num];

               loc_num++;
            }
         }

         detF = F[0]*F[3] - F[1]*F[2];
      }
      else //confflag == 2.  calculation of F^-1
      {
         double Finv[4] = {0.0, 0.0, 0.0, 0.0};

         loc_num = 0;
         for(jj=0; jj<=q; jj++)
         {
            for(ii=0; ii<=p; ii++)
            {
               EP = &(PP[ni+ii][nj+jj] );

               Finv[0] += EP->x * dN_dx[loc_num];
               Finv[2] += EP->x * dN_dy[loc_num];
               Finv[1] += EP->y * dN_dx[loc_num];
               Finv[3] += EP->y * dN_dy[loc_num];
               loc_num++;
            }
         }

         detF = 1.0 / (Finv[0]*Finv[3] - Finv[1]*Finv[2]);

         F[0] =   Finv[3] * detF ;
         F[1] = - Finv[1] * detF ;
         F[2] = - Finv[2] * detF ;
         F[3] =   Finv[0] * detF ;
      }


  return;
}
//


/*
void NurbsSURFACE::deformationGradient(int ni, int nj, bool confflag, double* dN_dx, double* dN_dy, double* F, double& detF)
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

      //cout << ni << '\t' << nj << endl;
      //cout << ngbf1 << '\t' << ngbf2 << endl;

      int  ii, jj, loc_num, temp;

      EPOINT *EP;

      F[0] = F[1] = F[2] = F[3] = 0.0;

      if(confflag) // calculation of F
      {
        int  ii, index;
        double  b1, b2;

         loc_num = 0;
         for(jj=0; jj<=q; jj++)
         {
            for(ii=0; ii<=p; ii++)
            {
               index = ngbf1*(nj+jj) + ni + ii;

               b1 = Values[0][index];
               b2 = Values[1][index];

               F[0] += b1 * dN_dx[loc_num];
               F[2] += b1 * dN_dy[loc_num];
               F[1] += b2 * dN_dx[loc_num];
               F[3] += b2 * dN_dy[loc_num];

               loc_num++;
            }
         }

         F[0] += 1.0;
         F[3] += 1.0;

         detF = F[0]*F[3] - F[1]*F[2];
      }
      else //confflag == 2.  calculation of F^-1
      {
         double Finv[4] = {0.0, 0.0, 0.0, 0.0};

         loc_num = 0;
         for(jj=0; jj<=q; jj++)
         {
            for(ii=0; ii<=p; ii++)
            {
               EP = &(PP[ni+ii][nj+jj] );

               Finv[0] += EP->x * dN_dx[loc_num];
               Finv[2] += EP->x * dN_dy[loc_num];
               Finv[1] += EP->y * dN_dx[loc_num];
               Finv[3] += EP->y * dN_dy[loc_num];
               loc_num++;
            }
         }

         detF = 1.0 / (Finv[0]*Finv[3] - Finv[1]*Finv[2]);

         F[0] =   Finv[3] * detF ;
         F[1] = - Finv[1] * detF ;
         F[2] = - Finv[2] * detF ;
         F[3] =   Finv[0] * detF ;
      }
      
  return;
}
*/



void NurbsSURFACE::deformationGradient2(int* tt, double* NN, double* dN_dx, double* dN_dy, MatrixXd& F, MatrixXd& b, double& pres)
{
    int  ii, index;
    double  b1, b2;

    F.setZero();
    for(ii=0;ii<nlbf;ii++)
    {
        index = tt[ii];

        b1 = Values[0][index];
        b2 = Values[1][index];

        F(0,0) += b1*dN_dx[ii];
        F(0,1) += b1*dN_dy[ii];
        F(1,0) += b2*dN_dx[ii];
        F(1,1) += b2*dN_dy[ii];
    }

    F(0,0) += 1.0;
    F(1,1) += 1.0;
    F(2,2) = 1.0;

    b = F*F.transpose();

    if(ndof > 2)
    {
       pres = 0.0;
       for(ii=0;ii<nlbf;ii++)
          pres += Values[2][tt[ii]] * NN[ii];
    }

  return;
}





void NurbsSURFACE::PlotValues(int  index)
{
    int ii, jj, n1, n2, ind1, ind2, ll, resln[2], nCol=9;
    resln[0] = resln[1] = 5;

    VectorArray<double>   uu, vv;

    if(p == 1)
      findunique(U, uu);
    else
      create_vector2(U, resln[0], uu);

    if(q == 1)
      findunique(V, vv);
    else
      create_vector2(V, resln[1], vv);

    cout << " AAAAAAAAAAA " << endl;

      vtkSmartPointer<vtkPoints>             points   =  vtkSmartPointer<vtkPoints>::New();
      vtkSmartPointer<vtkQuad>               quad     =  vtkSmartPointer<vtkQuad>::New();
      vtkSmartPointer<vtkFloatArray>        scalars   =   vtkSmartPointer<vtkFloatArray>::New();
      vtkSmartPointer<vtkScalarBarActor>    scalarBar =   vtkSmartPointer<vtkScalarBarActor>::New();
      vtkSmartPointer<vtkLookupTable>       hueLut    =   vtkSmartPointer<vtkLookupTable>::New();
      vtkSmartPointer<vtkDataSetMapper>     mapper1   =   vtkSmartPointer<vtkDataSetMapper>::New();
      vtkSmartPointer<vtkActor>             actor1    =   vtkSmartPointer<vtkActor>::New();


      vtkIdType pts[4], count;
      EPOINT  EP2;

      n1 = uu.n;    n2 = vv.n;
    
      scalars->SetNumberOfTuples(n1*n2);

      float   val;
      double  umn = 0.0, umx = 0.0;

      count=0; 

      for(jj=0;jj<n2;jj++)
      {
          for(ii=0;ii<n1;ii++)
          {
              EP2 = SurfacePoint(uu[ii], vv[jj]).CalcEuclid();

              points->InsertNextPoint(EP2.x, EP2.y, EP2.z);

              val = computeValue(1, uu[ii], vv[jj]) ;
              cout << " val " << val << endl;
              
              //val = EP2.z;

              if(val < umn) umn = val;
              if(val > umx) umx = val;
          
              scalars->SetTuple1(count++, val);
          }
      }


    for(jj=0;jj<n2-1;jj++)
    {
       ind1 = n1*jj;
       ind2 = n1*(jj+1);
       
       for(ii=0;ii<n1-1;ii++)
       {
          pts[0] = ind1+ii;          pts[1] = pts[0]+1;
          pts[3] = ind2+ii;          pts[2] = pts[3]+1;
          
          for(ll=0;ll<4;ll++)
             quad->GetPointIds()->SetId(ll, pts[ll]);
          
          plotvtk.uGrid->InsertNextCell(quad->GetCellType(), quad->GetPointIds());
       }
    }


    mapper1->ScalarVisibilityOn();
//    mapper1->SetScalarModeToUsePointData();
    mapper1->SetColorModeToMapScalars();
    
    mapper1->SetScalarRange(umn, umx);

    hueLut->SetHueRange(0.66667, 0.0);
    hueLut->SetTableRange(umn, umx);
//    hueLut->SetNumberOfColors(nCol);
//    hueLut->SetRampToLinear();
//    hueLut->SetScaleToLinear();
    hueLut->Build();
    
    vtkSmartPointer<vtkTextProperty>    txtprop       =   vtkSmartPointer<vtkTextProperty>::New();
    
    txtprop->SetColor(0.0,0.0,0.0);
    txtprop->BoldOn();
    txtprop->ItalicOff();


    mapper1->SetLookupTable(hueLut);
//    mapper1->UseLookupTableScalarRangeOn();

    scalarBar->SetLookupTable(hueLut);
//    scalars->SetLookupTable(hueLut);

    scalarBar->SetMaximumNumberOfColors(nCol);
    scalarBar->SetNumberOfLabels(nCol+1);
    scalarBar->SetTitleTextProperty(txtprop);
    scalarBar->SetLabelTextProperty(txtprop);
    scalarBar->SetWidth(0.12);
    scalarBar->SetLabelFormat("%10.4f");
    
    plotvtk.uGrid->SetPoints(points);

    plotvtk.uGrid->GetPointData()->SetScalars(scalars);

    mapper1->SetInputConnection(plotvtk.uGrid->GetProducerPort());
    //mapper1->SetInputData(plotvtk.uGrid);

    actor1->SetMapper(mapper1);

    plotvtk.rendr->AddActor(actor1);
    plotvtk.rendr->AddActor2D(scalarBar);
    plotvtk.rendr->ResetCamera();
    plotvtk.renWindow->Render();

  return;
}






void NurbsSURFACE::PlotElements(int col, bool PLOT_KNOT_LINES, int* resln)
{
      if(plotvtk.ActiveFlag)
      {
         PlotElementsVTK(col, PLOT_KNOT_LINES, resln);
         return;
      }

      VectorArray<double> uu, vv;

      if(p == 1)
        findunique(U, uu);
      else
        create_vector2(U, resln[0], uu);

      if(q == 1)
        findunique(V, vv);
      else
        create_vector2(V, resln[1], vv);

      int ii, jj;

      ListArray<ListArray<EPOINT> > S1;
      S1.setDim(uu.n);

      for(ii=0;ii<uu.n;ii++)
        S1[ii].setDim(vv.n);

      EPOINT *PP1, *PP2;

      // Calculate the Points on the surface
      for(ii=0;ii<uu.n;ii++)
      {
        for(jj=0;jj<vv.n;jj++)
        {
          S1[ii][jj] = SurfacePoint(uu[ii], vv[jj]).CalcEuclid();
        }
      }

      double x1d[3], x2d[3];

      // plot 'u' isolines
      for(ii=0;ii<uu.n;ii++)
      {
        if( IsKnot(&(U[0]), U.n ,uu[ii]) )
        {
          for(jj=0;jj<vv.n-1;jj++)
          {
            PP1 = &(S1[ii][jj]);
            PP2 = &(S1[ii][jj+1]);

            x1d[0] = PP1->x; x1d[1] = PP1->y; x1d[2] = PP1->z;
            x2d[0] = PP2->x; x2d[1] = PP2->y; x2d[2] = PP2->z;

            plot.line(x1d,x2d);
          }
        }
      }

      // plot 'v' isolines
      for(ii=0;ii<uu.n-1;ii++)
      {
        for(jj=0;jj<vv.n;jj++)
        {
          if( IsKnot(&(V[0]), V.n ,vv[jj]) )
          {
            PP1 = &(S1[ii][jj]);
            PP2 = &(S1[ii+1][jj]);

            x1d[0] = PP1->x; x1d[1] = PP1->y; x1d[2] = PP1->z;
            x2d[0] = PP2->x; x2d[1] = PP2->y; x2d[2] = PP2->z;

            plot.line(x1d,x2d);
          }
        }
      }

  return;

}







void NurbsSURFACE::PlotElementsVTK(int col, bool PLOT_KNOT_LINES, int* resln)
{
    int ii, jj, n1, n2, ind1, ind2, ll;
    double  color[3];

    VectorArray<double>   uu, vv;
    
    if(p == 1)
      findunique(U, uu);
    else
      create_vector2(U, resln[0], uu);

    if(q == 1)
      findunique(V, vv);
    else
      create_vector2(V, resln[1], vv);


    vtkSmartPointer<vtkPoints>             points   =  vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkUnstructuredGrid>   uGrid    =  vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkQuad>               quad     =  vtkSmartPointer<vtkQuad>::New();
    vtkSmartPointer<vtkDataSetMapper>      mapper1  =  vtkSmartPointer<vtkDataSetMapper>::New();
    vtkSmartPointer<vtkActor>              actor1   =  vtkSmartPointer<vtkActor>::New();

    vtkIdType pts[4];
    EPOINT  EP1;

    n1 = uu.n;    n2 = vv.n;
    
    //cout << " number of points " <<   plotvtk.uGrid->GetNumberOfPoints() << endl;

    for(jj=0;jj<n2;jj++)
    {
       for(ii=0;ii<n1;ii++)
       {
          EP1 = SurfacePoint(uu[ii], vv[jj]).CalcEuclid();
          points->InsertNextPoint(EP1.x, EP1.y, EP1.z);
       }
    }

    uGrid->SetPoints(points);

    for(jj=0;jj<n2-1;jj++)
    {
       ind1 = n1*jj;
       ind2 = n1*(jj+1);
       
       for(ii=0;ii<n1-1;ii++)
       {
          pts[0] = ind1+ii;          pts[1] = pts[0]+1;
          pts[3] = ind2+ii;          pts[2] = pts[3]+1;
          
          for(ll=0;ll<4;ll++)
             quad->GetPointIds()->SetId(ll, pts[ll]);
          
          uGrid->InsertNextCell(quad->GetCellType(), quad->GetPointIds());
       }
    }

    getColorValue(col, color);

    mapper1->SetInputConnection(uGrid->GetProducerPort());
    //mapper1->SetInputData(uGrid);
    
    actor1->SetMapper(mapper1);
    actor1->GetProperty()->SetColor(color[0], color[1], color[2]);
    actor1->GetProperty()->SetEdgeColor(0.0, 0.0, 0.0);
    actor1->GetProperty()->SetLineWidth(2.0);

    plotvtk.rendr->AddActor(actor1);
    
    if(PLOT_KNOT_LINES)
    {
        vtkSmartPointer<vtkUnstructuredGrid>  uGrid2     =  vtkSmartPointer<vtkUnstructuredGrid>::New();
        vtkSmartPointer<vtkLine>              line       =  vtkSmartPointer<vtkLine>::New();
        vtkSmartPointer<vtkDataSetMapper>     mapper2    =  vtkSmartPointer<vtkDataSetMapper>::New();
        vtkSmartPointer<vtkActor>             actor2     =  vtkSmartPointer<vtkActor>::New();

        uGrid2->SetPoints(points);
        
        vtkIdType  pt1, pt2;
        EPOINT  EP2;
        // constant 'u' knot lines
        for(ii=0;ii<n1;ii++)
        {
           if( IsKnot(&(U[0]), U.n ,uu[ii]) )
           {
              for(jj=0;jj<n2-1;jj++)
              {
                  pt1 = n1*jj+ii;
                  pt2 = n1*(jj+1)+ii;
                  
                  line->GetPointIds()->SetId(0,pt1);
                  line->GetPointIds()->SetId(1,pt2);
                  uGrid2->InsertNextCell(line->GetCellType(), line->GetPointIds());
              }
           }
        }
        // constant 'v' knot lines
        for(jj=0;jj<n2;jj++)
        {
           ind1 = n1*jj;
           if( IsKnot(&(V[0]), V.n ,vv[jj]) )
           {
              for(ii=0;ii<n1-1;ii++)
              {
                  pt1 = ind1+ii;
                  pt2 = pt1+1;
                  
                  line->GetPointIds()->SetId(0,pt1);
                  line->GetPointIds()->SetId(1,pt2);
                  uGrid2->InsertNextCell(line->GetCellType(), line->GetPointIds());
              }
           }
        }

        mapper2->SetInputConnection(uGrid2->GetProducerPort());
        //mapper2->SetInputData(uGrid2);

        actor2->SetMapper(mapper2);

        actor2->GetProperty()->SetEdgeColor(0.0, 0.0, 0.0); //(R,G,B)
    
        actor2->GetProperty()->SetLineWidth(2.0);
    
        actor2->GetProperty()->EdgeVisibilityOn();

        plotvtk.rendr->AddActor(actor2);
    } // if(PLOT_KNOT_LINES)

//    plotvtk.rendr->ResetCamera();

//    plotvtk.renWindow->Render();

  return;

}

