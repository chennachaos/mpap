
#include "GeomDataHBSplines.h"

#include "ComputerTime.h"
#include "TimeFunction.h"
#include "MpapTime.h"
#include "Functions.h"
#include "QuadratureUtil.h"
#include "ShapeFunctions.h"
#include "DistFunctions.h"
#include "Functions.h"
//#include "headersBoost.h"

#include "BasisFunctionsBSpline.h"

extern ComputerTime       computerTime;
extern MpapTime mpapTime;
extern List<TimeFunction> timeFunction;


GeomDataHBSplines::GeomDataHBSplines()
{
  degree[0] = degree[1] = degree[2] = 1;
  
  analyDBC = NULL;
}



GeomDataHBSplines::~GeomDataHBSplines()
{
  delete analyDBC;

  for(vector<DistanceFunction*>::iterator pObj = distFuncs.begin(); pObj != distFuncs.end(); ++pObj)
  {
    delete *pObj; // Note that this is deleting what pObj points to, which is a pointer
  }
  distFuncs.clear(); // Purge the contents so no one tries to delete them again

  /*
  for(int ii=0;ii<ROWS;ii++)
    {
      delete [] ders1[ii];
      delete [] ders2[ii];
      delete [] ders3[ii];
    }
    delete [] ders1;
    delete [] ders2;
    delete [] ders3;
  */
}



void GeomDataHBSplines::build()
{
    int size, ii, jj, kk, ind, ind1, ind2, ind3, lev, gp1, gp2, gp3;

    nlbf = totalNGP = 1;
    Jfull = 1.0;
    for(ii=0;ii<DIM;ii++)
    {
      ngbf[ii] = nelem[ii] + degree[ii];
      nGP[ii]  = FluidProps[0];
      dX[ii]   = gridLEN[ii]/nelem[ii];

      //as the parameter range is from 0.0 to 1.0, Jacobian in each direction is just the length of the domain in the respective direction
      Jacobian[ii] = gridLEN[ii];

      nlbf     *= (degree[ii]+1);
      totalNGP *= nGP[ii];
      Jfull    *= Jacobian[ii];
    }

    //cout << " Jacobian[0] " << Jacobian[0] << '\t' << Jacobian[1] << '\t' << Jacobian[2] << '\t'<< Jfull << endl;
    
    cout << " totalNGP " << nGP[0] << '\t' << nGP[1] << '\t' << totalNGP << endl;

    // compute subdivision coefficient matrices

    size = degree[0]+1;
        
    coeffLeft.resize(size,size);
    coeffRight.resize(size,size);

    GenerateCoeffMatrices(degree[0], coeffLeft, coeffRight);

    getGaussPoints1D(nGP[0], gausspoints1, gaussweights1);
    //getGaussPoints1D(nGP[0], gausspoints2, gaussweights2);
    //getGaussPoints1D(nGP[0], gausspoints3, gaussweights3);
    gausspoints2 = gausspoints1;
    gausspoints3 = gausspoints1;
    
    gaussweights2 = gaussweights1;
    gaussweights3 = gaussweights1;

    if(DIM == 1)
    {
      ind = nGP[0];
      gausspoints.resize(ind);
      gaussweights.resize(ind);
  
      ind=0;
      for(ii=0; ii<nGP[0]; ii++)
      {
        gausspoints[ind][0] = gausspoints1[ii];

        gaussweights[ind] = gaussweights1[ii];
        ind++;
      }
    }
    if(DIM == 2)
    {
      ind = nGP[0]*nGP[0];
      gausspoints.resize(ind);
      gaussweights.resize(ind);
  
      ind=0;
      for(jj=0; jj<nGP[0]; jj++)
      {
        for(ii=0; ii<nGP[0]; ii++)
        {
          gausspoints[ind][0] = gausspoints1[ii];
          gausspoints[ind][1] = gausspoints2[jj];

          gaussweights[ind] = gaussweights2[jj]*gaussweights1[ii];
          ind++;
        }
      }
    }
    if(DIM == 3)
    {
      ind = nGP[0]*nGP[0];
      gausspoints2D.resize(ind);
      gaussweights2D.resize(ind);
  
      ind=0;
      for(jj=0; jj<nGP[0]; jj++)
      {
        for(ii=0; ii<nGP[0]; ii++)
        {
          gausspoints2D[ind][0] = gausspoints1[ii];
          gausspoints2D[ind][1] = gausspoints2[jj];

          gaussweights2D[ind] = gaussweights2[jj]*gaussweights1[ii];
          ind++;
        }
      }

      ind = nGP[0]*nGP[0]*nGP[0];
      gausspoints.resize(ind);
      gaussweights.resize(ind);

      //cout << DIM << '\t' << ind << endl;
  
      ind=0;
      for(kk=0; kk<nGP[0]; kk++)
      {
        for(jj=0; jj<nGP[0]; jj++)
        {
          for(ii=0; ii<nGP[0]; ii++)
          {
            gausspoints[ind][0] = gausspoints1[ii];
            gausspoints[ind][1] = gausspoints2[jj];
            gausspoints[ind][2] = gausspoints3[kk];

            gaussweights[ind] = gaussweights3[kk]*gaussweights2[jj]*gaussweights1[ii];
            ind++;
          }
        }
      }
    }



    if(DIM > 1)
    {
        TensorProduct(coeffLeft,  coeffLeft,  coeffSW);
        TensorProduct(coeffLeft,  coeffRight, coeffSE);
        TensorProduct(coeffRight, coeffLeft,  coeffNW);
        TensorProduct(coeffRight, coeffRight, coeffNE);
        /*
        printMatrix(coeffLeft);
        printf("\n\n");
        printMatrix(coeffRight);
        printf("\n\n");
        printMatrix(coeffSW);
        printf("\n\n");
        printMatrix(coeffSE);
        printf("\n\n");
        printMatrix(coeffNW);
        printf("\n\n");
        printMatrix(coeffNE);
        printf("\n\n");
        */
    }
    if(DIM > 2)
    {
        TensorProduct(coeffLeft,  coeffSW, coeffSW_Back);
        TensorProduct(coeffLeft,  coeffSE, coeffSE_Back);
        TensorProduct(coeffLeft,  coeffNW, coeffNW_Back);
        TensorProduct(coeffLeft,  coeffNE, coeffNE_Back);

        TensorProduct(coeffRight, coeffSW, coeffSW_Front);
        TensorProduct(coeffRight, coeffSE, coeffSE_Front);
        TensorProduct(coeffRight, coeffNW, coeffNW_Front);
        TensorProduct(coeffRight, coeffNE, coeffNE_Front);
    }
    
    gpsLeft.push_back(-1.0);
    gpsRight.push_back(1.0);

    gwsLeft.push_back(1.0);
    gwsRight.push_back(1.0);

    /*
    printMatrix(coeffLeft);
    printf("\n\n");
    printMatrix(coeffRight);
    */
/*
    ROWS = 3;

    ders1 = new double*[ROWS];
    ders2 = new double*[ROWS];
    ders3 = new double*[ROWS];

    for(ii=0;ii<ROWS;ii++)
    {
       ders1[ii] = new double[degree[0]+1];
       ders2[ii] = new double[degree[1]+1];
       ders3[ii] = new double[degree[2]+1];
    }
*/

//
    double  du = 1.0/nelem[0];
    double  dv = 1.0/nelem[1];
    double  dw = 1.0/nelem[2];
    double  incr1, incr2, incr3, val1, val2, val3, val4, val5, val6, uu, vv, ww;

    shpfns.resize(8); // index 1 - level(s)
    
    myPoint  knotBegin, knotIncr, param;
    knotBegin.setZero();
    knotIncr = knotBegin;
    param = knotBegin;

    incr1 = du;    incr2 = dv;    incr3 = dw;

    for(lev=0;lev<8;lev++)
    {
      shpfns[lev].resize(totalNGP);

      val1  = 0.5*incr1;
      val2  = 0.5*incr1;

      val3  = 0.5*incr2;
      val4  = 0.5*incr2;

      val5  = 0.5*incr3;
      val6  = 0.5*incr3;
      
      knotIncr[0] = incr1;
      knotIncr[1] = incr2;
      knotIncr[2] = incr3;

      if(DIM == 2)
      {
        ind = 0;
        for(gp2=0;gp2<nGP[1];gp2++)
        {
          param[1] = val3 * gausspoints2[gp2] + val4;
        for(gp1=0;gp1<nGP[0];gp1++)
        {
          param[0] = val1 * gausspoints1[gp1] + val2;

          shpfns[lev][ind].Initialise(degree[0], degree[1]);

          computeBasisFunctions2D(knotBegin, knotIncr, param, shpfns[lev][ind].N, 
                                    shpfns[lev][ind].dN_dx, shpfns[lev][ind].dN_dy,
                                    shpfns[lev][ind].d2N_dx2, shpfns[lev][ind].d2N_dy2);

          //printVector(shpfns[lev][ind].N);
          //printf("\n\n\n");
          //printVector(shpfns[lev][ind].dN_dx);
          //printf("\n\n\n");
          //printVector(shpfns[lev][ind].dN_dy);
          //printf("\n\n\n");


          ind++;
        }
        }
      }
      if(DIM == 3)
      {
        ind = 0;
        for(gp3=0;gp3<nGP[2];gp3++)
        {
          param[2] = val5 * gausspoints2[gp3] + val6;
        for(gp2=0;gp2<nGP[1];gp2++)
        {
          param[1] = val3 * gausspoints2[gp2] + val4;
        for(gp1=0;gp1<nGP[0];gp1++)
        {
          param[0] = val1 * gausspoints1[gp1] + val2;

          shpfns[lev][ind].Initialise(degree[0], degree[1], degree[2]);

          computeBasisFunctions3D(knotBegin, knotIncr, param, shpfns[lev][ind].N, 
                        shpfns[lev][ind].dN_dx,   shpfns[lev][ind].dN_dy,   shpfns[lev][ind].dN_dz,
                        shpfns[lev][ind].d2N_dx2, shpfns[lev][ind].d2N_dy2, shpfns[lev][ind].d2N_dz2);
          ind++;
        }
        }
        }
      }

      incr1 *= 0.5;
      incr2 *= 0.5;
      incr3 *= 0.5;

      if(2 < 1)
      {
        printf("AAAAAAAAAAAAAAAAA \n");
        printVector(shpfns[lev][0].N);
        printf("\n\n");
        printVector(shpfns[lev][0].dN_dx);
        printf("\n\n");
        printVector(shpfns[lev][0].dN_dy);
        printf("\n\n");
        printVector(shpfns[lev][0].d2N_dx2);
        printf("\n\n");
        printVector(shpfns[lev][0].d2N_dy2);
        printf("\n\n");
      }
    }
    //printf("KKKKKKKKKKK \n");
//

    // set boundary normals
    // used to imposing boundary conditions on the boundaries of the background grid
    ////////////////////////////////////////

    boundaryNormals.resize(6);

    boundaryNormals[0][0] = -1.0; boundaryNormals[0][1] =  0.0; boundaryNormals[0][2] =  0.0;

    boundaryNormals[1][0] =  1.0; boundaryNormals[1][1] =  0.0; boundaryNormals[1][2] =  0.0;

    boundaryNormals[2][0] =  0.0; boundaryNormals[2][1] = -1.0; boundaryNormals[2][2] =  0.0;

    boundaryNormals[3][0] =  0.0; boundaryNormals[3][1] =  1.0; boundaryNormals[3][2] =  0.0;

    boundaryNormals[4][0] =  0.0; boundaryNormals[4][1] =  0.0; boundaryNormals[4][2] = -1.0;

    boundaryNormals[5][0] =  0.0; boundaryNormals[5][1] =  0.0; boundaryNormals[5][2] =  1.0;

    // set boundary Jacobians
    // used to imposing boundary conditions on the boundaries of the background grid
    ////////////////////////////////////////

    boundaryJacobians.resize(6);
    for(ii=0; ii<6; ii++)
      boundaryJacobians[ii].resize(10);


    if(DIM == 2)
    {
      val1 = 0.5 * dv * Jacobian[1] ;
      val2 = 0.5 * du * Jacobian[0] ;

      for(lev=0; lev<10; lev++)
      {
         boundaryJacobians[0][lev] = val1;
         boundaryJacobians[1][lev] = val1;

         boundaryJacobians[2][lev] = val2;
         boundaryJacobians[3][lev] = val2;

         val1 *= 0.5;
         val2 *= 0.5;
      }
    }

    if(DIM == 3)
    {
      val1 = 0.25 * dv * dw * Jacobian[1] * Jacobian[2];
      
      val2 = 0.25 * du * dw * Jacobian[0] * Jacobian[2];

      val3 = 0.25 * du * dv * Jacobian[0] * Jacobian[1];

      for(lev=0; lev<10; lev++)
      {
         boundaryJacobians[0][lev] = val1;
         boundaryJacobians[1][lev] = val1;

         boundaryJacobians[2][lev] = val2;
         boundaryJacobians[3][lev] = val2;

         boundaryJacobians[4][lev] = val3;
         boundaryJacobians[5][lev] = val3;

         val1 *= 0.25;
         val2 *= 0.25;
         val3 *= 0.25;
      }
    }

    // set boundary Gauss points and weights
    // used to imposing boundary conditions on the boundaries of the background grid
    ////////////////////////////////////////


    myPoint  ptTemp;
    double  wt;
    
    if(DIM == 2)
    {
       boundaryQuadrature2D.resize(4);

       ptTemp[2] = 0.0;

       //side 0:

              for(jj=0; jj<gausspoints2.size(); jj++)
              {
                ptTemp[1] = gausspoints2[jj];

                for(ii=0; ii<gpsLeft.size(); ii++)
                {
                  ptTemp[0] = gpsLeft[ii];

                  wt = gaussweights2[jj]*gwsLeft[ii];

                  boundaryQuadrature2D[0].gausspoints.push_back(ptTemp);
                  boundaryQuadrature2D[0].gaussweights.push_back(wt);
                }
              }

       //side 1:

              for(jj=0; jj<gausspoints2.size(); jj++)
              {
                ptTemp[1] = gausspoints2[jj];

                for(ii=0; ii<gpsRight.size(); ii++)
                {
                  ptTemp[0] = gpsRight[ii];

                  wt = gaussweights2[jj]*gwsLeft[ii];

                  boundaryQuadrature2D[1].gausspoints.push_back(ptTemp);
                  boundaryQuadrature2D[1].gaussweights.push_back(wt);
                }
              }

       //side 2:

              for(jj=0; jj<gpsLeft.size(); jj++)
              {
                ptTemp[1] = gpsLeft[jj];

                for(ii=0; ii<gausspoints1.size(); ii++)
                {
                  ptTemp[0] = gausspoints1[ii];

                  wt = gwsLeft[jj]*gaussweights1[ii];

                  boundaryQuadrature2D[2].gausspoints.push_back(ptTemp);
                  boundaryQuadrature2D[2].gaussweights.push_back(wt);
                }
              }

       //side 3:

              for(jj=0; jj<gpsRight.size(); jj++)
              {
                ptTemp[1] = gpsRight[jj];

                for(ii=0; ii<gausspoints1.size(); ii++)
                {
                  ptTemp[0] = gausspoints1[ii];

                  wt = gwsRight[jj]*gaussweights1[ii];

                  boundaryQuadrature2D[3].gausspoints.push_back(ptTemp);
                  boundaryQuadrature2D[3].gaussweights.push_back(wt);
                }
              }
    }


    if(DIM == 3)
    {
       boundaryQuadrature3D.resize(6);

       //side 0:

            for(kk=0; kk<gausspoints3.size(); kk++)
            {
              ptTemp[2] = gausspoints3[kk];

              for(jj=0; jj<gausspoints2.size(); jj++)
              {
                ptTemp[1] = gausspoints2[jj];

                for(ii=0; ii<gpsLeft.size(); ii++)
                {
                  ptTemp[0] = gpsLeft[ii];

                  wt = gaussweights3[kk]*gaussweights2[jj]*gwsLeft[ii];
                  
                  boundaryQuadrature3D[0].gausspoints.push_back(ptTemp);
                  boundaryQuadrature3D[0].gaussweights.push_back(wt);
                }
              }
            }

       //side 1:

            for(kk=0; kk<gausspoints3.size(); kk++)
            {
              ptTemp[2] = gausspoints3[kk];

              for(jj=0; jj<gausspoints2.size(); jj++)
              {
                ptTemp[1] = gausspoints2[jj];

                for(ii=0; ii<gpsRight.size(); ii++)
                {
                  ptTemp[0] = gpsRight[ii];

                  wt = gaussweights3[kk]*gaussweights2[jj]*gwsRight[ii];

                  boundaryQuadrature3D[1].gausspoints.push_back(ptTemp);
                  boundaryQuadrature3D[1].gaussweights.push_back(wt);
                }
              }
            }

       //side 2:

            for(kk=0; kk<gausspoints3.size(); kk++)
            {
              ptTemp[2] = gausspoints3[kk];

              for(jj=0; jj<gpsLeft.size(); jj++)
              {
                ptTemp[1] = gpsLeft[jj];

                for(ii=0; ii<gausspoints1.size(); ii++)
                {
                  ptTemp[0] = gausspoints1[ii];

                  wt = gaussweights3[kk]*gwsLeft[jj]*gaussweights1[ii];

                  boundaryQuadrature3D[2].gausspoints.push_back(ptTemp);
                  boundaryQuadrature3D[2].gaussweights.push_back(wt);
                }
              }
            }

       //side 3:

            for(kk=0; kk<gausspoints3.size(); kk++)
            {
              ptTemp[2] = gausspoints3[kk];

              for(jj=0; jj<gpsRight.size(); jj++)
              {
                ptTemp[1] = gpsRight[jj];

                for(ii=0; ii<gausspoints1.size(); ii++)
                {
                  ptTemp[0] = gausspoints1[ii];

                  wt = gaussweights3[kk]*gwsRight[jj]*gaussweights1[ii];

                  boundaryQuadrature3D[3].gausspoints.push_back(ptTemp);
                  boundaryQuadrature3D[3].gaussweights.push_back(wt);
                }
              }
            }


       //side 4:

            for(kk=0; kk<gpsLeft.size(); kk++)
            {
              ptTemp[2] = gpsLeft[kk];

              for(jj=0; jj<gausspoints2.size(); jj++)
              {
                ptTemp[1] = gausspoints2[jj];

                for(ii=0; ii<gausspoints1.size(); ii++)
                {
                  ptTemp[0] = gausspoints1[ii];

                  wt = gwsLeft[kk]*gaussweights2[jj]*gaussweights1[ii];

                  boundaryQuadrature3D[4].gausspoints.push_back(ptTemp);
                  boundaryQuadrature3D[4].gaussweights.push_back(wt);
                }
              }
            }

       //side 5:

            for(kk=0; kk<gpsRight.size(); kk++)
            {
              ptTemp[2] = gpsRight[kk];

              for(jj=0; jj<gausspoints2.size(); jj++)
              {
                ptTemp[1] = gausspoints2[jj];

                for(ii=0; ii<gausspoints1.size(); ii++)
                {
                  ptTemp[0] = gausspoints1[ii];

                  wt = gwsRight[kk]*gaussweights2[jj]*gaussweights1[ii];

                  boundaryQuadrature3D[5].gausspoints.push_back(ptTemp);
                  boundaryQuadrature3D[5].gaussweights.push_back(wt);
                }
              }
            }
    }

  ////////////////////////////////////////
  // analytical function for the Dirichlet boundary conditions
  //
  ////////////////////////////////////////
  analyDBC = (Function*) new PoissonEx3;
  //analyDBC = (Function*) new PoissonInterfaceEx4;


  //PoissonEx1 analy;
  //PoissonEx3 analy;
  //PoissonInterfaceEx3  analy(1.0, 6.0, 0.0);

  ////////////////////////////////////////
  // generate polygons for the immersed boundaries
  //
  ////////////////////////////////////////
  
  DistanceFunction  *distF1 = new Circle(0.0, 0.0, 0.5);
  Circle distF(0.0, 0.0, 0.5);
  
  distFuncs.push_back(distF1);

  return;
}




void GeomDataHBSplines::initialise(int size1, int size2, int size3, int size4)
{
  return;
}



void GeomDataHBSplines::printSelf()
{
//   cout << " Degree and Jacobian " << degree[0] << '\t' << Jfull << endl;
   return;
}


void GeomDataHBSplines::computeBasisFunctions1D(const myPoint& start, const myPoint& incr, const myPoint& param, VectorXd& N)
{
  HB_BasisFuns(degree[0], start[0], incr[0], param[0], &N(0));
  
  return;
}



void GeomDataHBSplines::computeBasisFunctions1D(const myPoint& start, const myPoint& incr, const myPoint& param, VectorXd& N, VectorXd& dN_dx)
{
    int  ROWS = 2, ii;

    double** ders = new double*[ROWS];

    for(ii=0;ii<ROWS;ii++)
       ders[ii] = new double[degree[0]+1];

    HB_DersBasisFuns(degree[0], start[0], incr[0], param[0], 1, ders);

    //  dx_du = J;
    //  du_dx = 1.0/dx_du

    double  du_dx = 1.0/Jacobian[0];
    
    for(ii=0; ii<=degree[0]; ii++)
    {
        N[ii]     =  ders[0][ii];
        dN_dx[ii] =  ders[1][ii] * du_dx;
    }

    for(ii=0;ii<ROWS;ii++)
      delete [] ders[ii];

    delete [] ders;
    
    return;
}



void GeomDataHBSplines::computeBasisFunctions1D(const myPoint& start, const myPoint& incr, const myPoint& param, VectorXd& N, VectorXd& dN_dx, VectorXd& d2N_dx2)
{
    int ROWS = 3, ii, Un;

    double** ders = new double*[ROWS];

    for(ii=0;ii<ROWS;ii++)
       ders[ii] = new double[degree[0]+1];

    HB_DersBasisFuns(degree[0], start[0], incr[0], param[0], 2, ders);

    //  dx_du = J; and     d2x_du2 = 0.0;
    //  du_dx = 1.0/dx_du
    //  d2u_dx2 = -d2x_du2/pow((dx_du),3) = 0.0
    // So,

    double  du_dx = 1.0/Jacobian[0];
    double  temp  = du_dx * du_dx;

    for(ii=0; ii<=degree[0]; ii++)
    {
        N[ii]       = ders[0][ii];
        dN_dx[ii]   = ders[1][ii] * du_dx;
        d2N_dx2[ii] = ders[2][ii] * temp;
    }

    for(ii=0;ii<ROWS;ii++)
      delete [] ders[ii];

    delete [] ders;

    return;
}



void GeomDataHBSplines::computeBasisFunctions1D(const myPoint& start, const myPoint& incr, const myPoint& param, VectorXd& N, VectorXd& dN_dx, VectorXd& d2N_dx2, VectorXd& d3N_dx3)
{
    int ROWS = 4, ii, Un;

    double** ders = new double*[ROWS];

    for(ii=0;ii<ROWS;ii++)
       ders[ii] = new double[degree[0]+1];

    HB_DersBasisFuns(degree[0], start[0], incr[0], param[0], 3, ders);

    //  dx_du = J; and     d2x_du2 = 0.0;
    //  du_dx = 1.0/dx_du
    //  d2u_dx2 = -d2x_du2/pow((dx_du),3) = 0.0
    // So,

    double  du_dx = 1.0/Jacobian[0];
    double  temp  = du_dx * du_dx;
    double  temp2 = temp * du_dx;

    for(ii=0; ii<=degree[0]; ii++)
    {
        N[ii]       = ders[0][ii];
        dN_dx[ii]   = ders[1][ii] * du_dx;
        d2N_dx2[ii] = ders[2][ii] * temp;
        d3N_dx3[ii] = ders[3][ii] * temp2;
    }

    for(ii=0;ii<ROWS;ii++)
      delete [] ders[ii];

    delete [] ders;

    return;
}




void GeomDataHBSplines::computeBasisFunctions2D(const myPoint& start, const myPoint& incr, const myPoint& param, VectorXd& N)
{
    int  ii, jj, count;

    vector<double>  N1(degree[0]+1), N2(degree[1]+1);

    HB_BasisFuns(degree[0], start[0], incr[0], param[0], &N1[0]);
    HB_BasisFuns(degree[1], start[1], incr[1], param[1], &N2[0]);

    count = 0;
    for(jj=0; jj<=degree[1]; jj++)
    {
       for(ii=0; ii<=degree[0]; ii++)
       {
          N[count++]  =  N1[ii] * N2[jj];
       }
    }

    return;
}



void GeomDataHBSplines::computeBasisFunctions2D(const myPoint& start, const myPoint& incr, const myPoint& param, 
					   VectorXd& N, VectorXd& dN_dx, VectorXd& dN_dy)
{
    //cout << " degree[0] " << degree[0] << '\t' << degree[1] << endl;
    int ROWS = 2, ii, jj, count;

    double** ders1 = new double*[ROWS];
    double** ders2 = new double*[ROWS];

    for(ii=0;ii<ROWS;ii++)
    {
       ders1[ii] = new double[degree[0]+1];
       ders2[ii] = new double[degree[1]+1];
    }
    //cout << " degree[0] " << degree[0] << '\t' << degree[1] << endl;
    HB_DersBasisFuns(degree[0], start[0], incr[0], param[0], 1, ders1);
    HB_DersBasisFuns(degree[1], start[1], incr[1], param[1], 1, ders2);
    //cout << " 222222222222222 " << endl;
    double  fact3, fact4;
    
    count = 0;
    for(jj=0; jj<=degree[1]; jj++)
    {
      fact3 = ders2[0][jj]/Jacobian[0];
      fact4 = ders2[1][jj]/Jacobian[1];

      for(ii=0; ii<=degree[0]; ii++)
      {
          N[count]       =  ders1[0][ii] * ders2[0][jj];

          dN_dx[count]   =  ders1[1][ii] * fact3;
          dN_dy[count]   =  ders1[0][ii] * fact4;
          count++;
      }
    }
    //cout << " 222222222222222 " << endl;
    for(ii=0;ii<ROWS;ii++)
    {
      delete [] ders1[ii];
      delete [] ders2[ii];
    }
    delete [] ders1;
    delete [] ders2;

    return;
}




void GeomDataHBSplines::computeBasisFunctionsGhostPenalty2D(const myPoint& start, const myPoint& incr, const myPoint& param, 
					   VectorXd& N, VectorXd& dN_dx, VectorXd& dN_dy)
{
  // derivatives of the basis functions
  // for ghost-penalty terms
  // 
  // compute the highest non-zero derivatives 
  // for the corresponding degree of B-Splines
  // 

    int ROWS = (degree[0]+1), ii, jj, count;

    double** ders1 = new double*[ROWS];
    double** ders2 = new double*[ROWS];

    for(ii=0;ii<ROWS;ii++)
    {
       ders1[ii] = new double[degree[0]+1];
       ders2[ii] = new double[degree[1]+1];
    }

    HB_DersBasisFuns(degree[0], start[0], incr[0], param[0], degree[0], ders1);
    HB_DersBasisFuns(degree[1], start[1], incr[1], param[1], degree[0], ders2);

    double  fact1, fact2;

    count = 0;
    for(jj=0; jj<=degree[1]; jj++)
    {
      fact1 = ders2[0][jj];
      fact2 = ders2[degree[0]][jj];

      for(ii=0; ii<=degree[0]; ii++)
      {
        dN_dx[count] = fact1 * ders1[degree[0]][ii] ;
        dN_dy[count] = fact2 * ders1[0][ii] ;

        count++;
      }
    }

    fact1 = fact2 = 1.0;
    for(ii=0; ii<degree[0]; ii++)
    {
      fact1 /= Jacobian[0];
      fact2 /= Jacobian[1];
    }

    count = (degree[0]+1)*(degree[1]+1);

    for(ii=0; ii<count; ii++)
    {
      dN_dx[ii] *= fact1 ;
      dN_dy[ii] *= fact2 ;
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



void GeomDataHBSplines::computeBasisFunctionsGhostPenalty3D(const myPoint& start, const myPoint& incr, const myPoint& param, 
					   VectorXd& N, VectorXd& dN_dx, VectorXd& dN_dy, VectorXd& dN_dz)
{
  // derivatives of the basis functions
  // for ghost-penalty terms
  // 
  // compute the highest non-zero derivatives 
  // for the corresponding degree of B-Splines
  // 

    int ROWS = (degree[0]+1), ii, jj, kk, count;

    double** ders1 = new double*[ROWS];
    double** ders2 = new double*[ROWS];
    double** ders3 = new double*[ROWS];

    for(ii=0;ii<ROWS;ii++)
    {
       ders1[ii] = new double[degree[0]+1];
       ders2[ii] = new double[degree[1]+1];
       ders3[ii] = new double[degree[2]+1];
    }

    HB_DersBasisFuns(degree[0], start[0], incr[0], param[0], degree[0], ders1);
    HB_DersBasisFuns(degree[1], start[1], incr[1], param[1], degree[1], ders2);
    HB_DersBasisFuns(degree[2], start[2], incr[2], param[2], degree[2], ders3);

    double  a1, a2, a3, b1, b2, b3;

    count = 0;
    for(kk=0;kk<=degree[2];kk++)
    {
      a1 = ders3[0][kk];
      a2 = ders3[0][kk];
      a3 = ders3[degree[0]][kk];

      for(jj=0;jj<=degree[1];jj++)
      {
        b1 = a1 * ders2[0][jj];
        b2 = a2 * ders2[degree[0]][jj];
        b3 = a3 * ders2[0][jj];

        for(ii=0;ii<=degree[0];ii++)
        {
          dN_dx[count]   = b1 * ders1[degree[0]][ii];
          dN_dy[count]   = b2 * ders1[0][ii];
          dN_dz[count]   = b3 * ders1[0][ii];

          count++;
        }
      }
    }

    a1 = a2 = a3 = 1.0;
    for(ii=0; ii<degree[0]; ii++)
    {
      a1 /= Jacobian[0];
      a2 /= Jacobian[1];
      a3 /= Jacobian[2];
    }

    count = (degree[0]+1)*(degree[1]+1)*(degree[2]+1);

    for(ii=0; ii<count; ii++)
    {
      dN_dx[ii] *= a1 ;
      dN_dy[ii] *= a2 ;
      dN_dz[ii] *= a3 ;
    }

    for(ii=0;ii<ROWS;ii++)
    {
      delete [] ders1[ii];
      delete [] ders2[ii];
      delete [] ders3[ii];
    }
    delete [] ders1;
    delete [] ders2;
    delete [] ders3;

    return;
}





void GeomDataHBSplines::computeBasisFunctions2D(const myPoint& start, const myPoint& incr, const myPoint& param, 
                              VectorXd& N, VectorXd& dN_dx, VectorXd& dN_dy, VectorXd& d2N_dx2, VectorXd& d2N_dy2)
{
    int ROWS = 3, ii, jj, count;

    double** ders1 = new double*[ROWS];
    double** ders2 = new double*[ROWS];

    for(ii=0;ii<ROWS;ii++)
    {
       ders1[ii] = new double[degree[0]+1];
       ders2[ii] = new double[degree[1]+1];
    }

    HB_DersBasisFuns(degree[0], start[0], incr[0], param[0], 2, ders1);
    HB_DersBasisFuns(degree[1], start[1], incr[1], param[1], 2, ders2);

    double  fact1, fact2, fact3, fact4, fact5, fact6;
    
    fact1 = Jacobian[0]*Jacobian[0];
    fact2 = Jacobian[1]*Jacobian[1];

    count = 0;
    for(jj=0; jj<=degree[1]; jj++)
    {
       fact3 = ders2[0][jj]/Jacobian[0];
       fact4 = ders2[1][jj]/Jacobian[1];

       fact5 = ders2[0][jj]/fact1;
       fact6 = ders2[2][jj]/fact2;
       
       for(ii=0; ii<=degree[0]; ii++)
       {
          N[count]       =  ders1[0][ii] * ders2[0][jj];
       
          dN_dx[count]   =  ders1[1][ii] * fact3;
          dN_dy[count]   =  ders1[0][ii] * fact4;

          d2N_dx2[count] =  ders1[2][ii] * fact5;
          d2N_dy2[count] =  ders1[0][ii] * fact6;

          count++;
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





void GeomDataHBSplines::computeBasisFunctions2D(const myPoint& start, const myPoint& incr, const myPoint& param, 
                              VectorXd& d3N_dx3, VectorXd& d3N_dy3, VectorXd& d3N_dxdy2, VectorXd& d3N_dx2dy)
{
    int ROWS = 4, ii, jj, count;

    double** ders1 = new double*[ROWS];
    double** ders2 = new double*[ROWS];

    for(ii=0;ii<ROWS;ii++)
    {
       ders1[ii] = new double[degree[0]+1];
       ders2[ii] = new double[degree[1]+1];
    }

    HB_DersBasisFuns(degree[0], start[0], incr[0], param[0], 3, ders1);
    HB_DersBasisFuns(degree[1], start[1], incr[1], param[1], 3, ders2);

    double  fact1, fact2, fact3, fact4, fact5, fact6;
    double  Jx2, Jy2, JxJy, Jx2Jy, JxJy2, Jx3, Jy3;
    
    Jx3 = Jacobian[0]*Jacobian[0]*Jacobian[0];
    Jy3 = Jacobian[1]*Jacobian[1]*Jacobian[1];
    
    Jx2Jy = Jacobian[0]*Jacobian[0]*Jacobian[1];
    JxJy2 = Jacobian[0]*Jacobian[1]*Jacobian[1];

    count = 0;
    for(jj=0; jj<=degree[1]; jj++)
    {
       for(ii=0; ii<=degree[0]; ii++)
       {
          d3N_dx3[count]    =  ders2[0][jj] * ders1[3][ii] / Jx3;
          d3N_dy3[count]    =  ders2[3][jj] * ders1[0][ii] / Jy3;

          d3N_dxdy2[count]  =  ders2[2][jj] * ders1[1][ii] / JxJy2;
          d3N_dx2dy[count]  =  ders2[1][jj] * ders1[2][ii] / Jx2Jy;

          count++;
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




void  GeomDataHBSplines::computeBasisFunctions3D(const myPoint& start, const myPoint& incr, const myPoint& param, VectorXd& N)
{
    int  ii, jj, kk, count;

    vector<double>  N1(degree[0]+1), N2(degree[1]+1), N3(degree[2]+1);

    HB_BasisFuns(degree[0], start[0], incr[0], param[0], &N1[0]);
    HB_BasisFuns(degree[1], start[1], incr[1], param[1], &N2[0]);
    HB_BasisFuns(degree[2], start[2], incr[2], param[2], &N3[0]);

    count = 0;
    for(kk=0;kk<=degree[2];kk++)
    {
      for(jj=0;jj<=degree[1];jj++)
      {
         for(ii=0;ii<=degree[0];ii++)
         {
           N[count++]  =  N1[ii] * N2[jj] * N3[kk];
         }
      }
    }

  return;
}

void  GeomDataHBSplines::computeBasisFunctions3D(const myPoint& start, const myPoint& incr, const myPoint& param, VectorXd& N, VectorXd& dN_dx, VectorXd& dN_dy, VectorXd& dN_dz)
{
    int ROWS = 2, ii, jj, kk, count;

    double** ders1 = new double*[ROWS];
    double** ders2 = new double*[ROWS];
    double** ders3 = new double*[ROWS];

    for(ii=0;ii<ROWS;ii++)
    {
       ders1[ii] = new double[degree[0]+1];
       ders2[ii] = new double[degree[1]+1];
       ders3[ii] = new double[degree[2]+1];
    }

    HB_DersBasisFuns(degree[0], start[0], incr[0], param[0], 1, ders1);
    HB_DersBasisFuns(degree[1], start[1], incr[1], param[1], 1, ders2);
    HB_DersBasisFuns(degree[2], start[2], incr[2], param[2], 1, ders3);

    double  a1, a2, a3, b1, b2, b3, b4;
    
    count = 0;
    for(kk=0;kk<=degree[2];kk++)
    {
      a1 = ders3[0][kk]/Jacobian[0];
      a2 = ders3[0][kk]/Jacobian[1];
      a3 = ders3[1][kk]/Jacobian[2];

      for(jj=0;jj<=degree[1];jj++)
      {
        b1 = ders2[0][jj]*a1;
        b2 = ders2[1][jj]*a2;
        b3 = ders2[0][jj]*a3;
        b4 = ders2[0][jj]*ders3[0][kk];

        for(ii=0;ii<=degree[0];ii++)
        {
          N[count]       =  ders1[0][ii] * b4;

          dN_dx[count]   =  ders1[1][ii] * b1;
          dN_dy[count]   =  ders1[0][ii] * b2;
          dN_dz[count]   =  ders1[0][ii] * b3;
          count++;
        }
      }
    }
    
    for(ii=0;ii<ROWS;ii++)
    {
      delete [] ders1[ii];
      delete [] ders2[ii];
      delete [] ders3[ii];
    }
    delete [] ders1;
    delete [] ders2;
    delete [] ders3;

  return;
}



void  GeomDataHBSplines::computeBasisFunctions3D(const myPoint& start, const myPoint& incr, const myPoint& param, VectorXd& N, VectorXd& dN_dx, VectorXd& dN_dy, VectorXd& dN_dz,
		VectorXd& d2N_dx2, VectorXd& d2N_dy2, VectorXd& d2N_dz2)
{
    int ROWS = 3, ii, jj, kk, count;

    double** ders1 = new double*[ROWS];
    double** ders2 = new double*[ROWS];
    double** ders3 = new double*[ROWS];

    for(ii=0;ii<ROWS;ii++)
    {
       ders1[ii] = new double[degree[0]+1];
       ders2[ii] = new double[degree[1]+1];
       ders3[ii] = new double[degree[2]+1];
    }

    HB_DersBasisFuns(degree[0], start[0], incr[0], param[0], 2, ders1);
    HB_DersBasisFuns(degree[1], start[1], incr[1], param[1], 2, ders2);
    HB_DersBasisFuns(degree[2], start[2], incr[2], param[2], 2, ders3);

    double  a1, a2, a3, a4, a5, a6, a7, b1, b2, b3, b4, b5, b6, b7, c1, c2, c3;

    c1 = Jacobian[0]*Jacobian[0];
    c2 = Jacobian[1]*Jacobian[1];
    c3 = Jacobian[2]*Jacobian[2];

    count = 0;
    for(kk=0;kk<=degree[2];kk++)
    {
      a1 = ders3[0][kk]/Jacobian[0];
      a2 = ders3[0][kk]/Jacobian[1];
      a3 = ders3[1][kk]/Jacobian[2];

      a5 = ders3[0][kk]/c1;
      a6 = ders3[0][kk]/c2;
      a7 = ders3[2][kk]/c3;

      for(jj=0;jj<=degree[1];jj++)
      {
        b1 = a1 * ders2[0][jj];
        b2 = a2 * ders2[1][jj];
        b3 = a3 * ders2[0][jj];
        b4 = ders3[0][kk] * ders2[0][jj];

        b5 = a5 * ders2[0][jj];
        b6 = a6 * ders2[2][jj];
        b7 = a7 * ders2[0][jj];

        for(ii=0;ii<=degree[0];ii++)
        {
          N[count]       = b4 * ders1[0][ii];

          dN_dx[count]   = b1 * ders1[1][ii];
          dN_dy[count]   = b2 * ders1[0][ii];
          dN_dz[count]   = b3 * ders1[0][ii];

          d2N_dx2[count] = b5 * ders1[2][ii];
          d2N_dy2[count] = b6 * ders1[0][ii];
          d2N_dz2[count] = b7 * ders1[0][ii];

          count++;
        }
      }
    }
    
    for(ii=0;ii<ROWS;ii++)
    {
      delete [] ders1[ii];
      delete [] ders2[ii];
      delete [] ders3[ii];
    }
    delete [] ders1;
    delete [] ders2;
    delete [] ders3;

  return;
}




void  GeomDataHBSplines::reset()
{
  return;
}



int  GeomDataHBSplines::doIntersect2D(AABB& bbTemp, bool flag, vector<int>& vecTemp, vector<myPoint>& ptOut, vector<int>& domVec)
{
  domVec.clear();
  int domTemp = 0;
  bool  ff=false;
  for(int bb=0; bb<immSolidPtrs.size(); bb++)
  {
    if( immSolidPtrs[bb]->doAABBintersect(bbTemp) )
    {
      domTemp = immSolidPtrs[bb]->doIntersect2D(bbTemp, flag, vecTemp, ptOut) ;

      if(domTemp == -1)
      {
        domVec.push_back(0);
        domVec.push_back(bb+1);

        ff = true;
      }
      else if(domTemp != 0)
      {
        domVec.push_back(domTemp);
        ff = false;
      }
      
      //cout << " ff = " << ff << '\t' << domTemp << endl;
      
      if(ff)
        break;
    }
  }

  if(domVec.empty())
    domVec.push_back(0);

  findUnique(domVec);

  return  1;
}



int  GeomDataHBSplines::doIntersect2Dfor3D(int sideTemp, double coord3, AABB& bbTemp, bool flag, vector<int>& vecTemp, vector<myPoint>& ptOut, vector<int>& domVec)
{
  domVec.clear();
  int domTemp = 0;
  bool  ff=false;
  for(int bb=0; bb<immSolidPtrs.size(); bb++)
  {
    //if( immSolidPtrs[bb]->doAABBintersect(bbTemp) )
    //{
      domTemp = immSolidPtrs[bb]->doIntersect2Dfor3D(sideTemp, coord3, bbTemp, flag, vecTemp, ptOut) ;

      if(domTemp == -1)
      {
        domVec.push_back(0);
        domVec.push_back(bb+1);

        ff = true;
      }
      else if(domTemp != 0)
      {
        domVec.push_back(domTemp);
        ff = false;
      }
      
      //cout << " ff = " << ff << '\t' << domTemp << endl;
      
      if(ff)
        break;
    //}
  }

  if(domVec.empty())
    domVec.push_back(0);

  findUnique(domVec);

  return  1;
}



int  GeomDataHBSplines::doIntersect3D(AABB& bbTemp, bool flag, vector<int>& vecTemp, vector<myPoint>& ptOut, vector<int>& domVec)
{
  domVec.clear();
  int domTemp = 0;
  bool  ff=false;

  for(int bb=0; bb<immSolidPtrs.size(); bb++)
  {
    if( immSolidPtrs[bb]->doAABBintersect(bbTemp) )
    {
      domTemp = immSolidPtrs[bb]->doIntersect3D(bbTemp, flag, vecTemp, ptOut) ;

      if(domTemp == -1)
      {
        domVec.push_back(0);
        domVec.push_back(bb+1);

        ff = true;
      }
      else if(domTemp != 0)
      {
        domVec.push_back(domTemp);
        ff = false;
      }
      
      //cout << " ff = " << ff << '\t' << domTemp << endl;
      
      if(ff)
        break;
    }
  }

  if(domVec.empty())
    domVec.push_back(0);

  findUnique(domVec);

  return  1;
}




int  GeomDataHBSplines::within(myPoint& ptTemp)
{
  int domainNum = 0;
  for(int bb=0; bb<immSolidPtrs.size(); bb++)
  {
    if( immSolidPtrs[bb]->within(ptTemp) )
    {
      domainNum = (bb+1) ;
      break;
    }
  }

  return  domainNum;
}








