
#include "GeomDataHBSplines.h"
#include "ComputerTime.h"
#include "TimeFunction.h"
#include "MpapTime.h"
#include "Functions.h"
#include "QuadratureUtil.h"
#include "ShapeFunctions.h"
#include "DistFunctions.h"
#include "Functions.h"
#include "BasisFunctionsBSpline.h"



void GeomDataHBSplines::getBoundaryNormal1D(int side, myPoint& normal)
{
    switch(side)
    {
        case 0:
          normal[0] = -1.0;
        break;

        case 1:
          normal[0] =  1.0;
        break;

        default :

             cout << " Invalid 'side' value in GeomDataHBSplines::getBoundaryNormal1D " << endl;
        break;
    } //switch(side)

    return;
}

void GeomDataHBSplines::getBoundaryNormal2D(int side, myPoint& normal)
{
    switch(side)
    {
        case 0:
          normal[0] = -1.0; normal[1] =  0.0;
        break;

        case 1:
          normal[0] =  1.0; normal[1] =  0.0;
        break;

        case 2:
          normal[0] =  0.0; normal[1] = -1.0;
        break;

        case 3:
          normal[0] =  0.0; normal[1] =  1.0;
        break;

        default :

            cout << " Invalid 'side' value in GeomDataHBSplines::getBoundaryNormal2D " << endl;
        break;
    } //switch(side)

    return;
}


void GeomDataHBSplines::getBoundaryNormal3D(int side, myPoint& normal)
{
    switch(side)
    {
        case 0:
          normal[0] = -1.0; normal[1] =  0.0; normal[2] =  0.0;
        break;

        case 1:
          normal[0] =  1.0; normal[1] =  0.0; normal[2] =  0.0;
        break;

        case 2:
          normal[0] =  0.0; normal[1] = -1.0; normal[2] =  0.0;
        break;

        case 3:
          normal[0] =  0.0; normal[1] =  1.0; normal[2] =  0.0;
        break;

        case 4:
          normal[0] =  0.0; normal[1] =  0.0; normal[2] = -1.0;
        break;

        case 5:
          normal[0] =  0.0; normal[1] =  0.0; normal[2] =  1.0;
        break;

        default :
            cout << " Invalid 'side' value in GeomDataHBSplines::getBoundaryNormal3D " << endl;
        break;
    } //switch(side)

    return;
}


void  GeomDataHBSplines::setBoundaryGPs1D(int side, vector<double>& boundaryGPs1, vector<double>& boundaryGWs1 )
{
    return;
}



void  GeomDataHBSplines::setBoundaryGPs2D(int side, vector<double>& boundaryGPs1, vector<double>& boundaryGWs1, vector<double>& boundaryGPs2, vector<double>& boundaryGWs2)
{
    switch(side)
    {
        case 0:

            boundaryGPs1 = gpsLeft;
            boundaryGWs1 = gwsLeft;

            boundaryGPs2 = gausspoints2;
            boundaryGWs2 = gaussweights2;

        break;

        case 1:

            boundaryGPs1 = gpsRight;
            boundaryGWs1 = gwsRight;

            boundaryGPs2 = gausspoints2;
            boundaryGWs2 = gaussweights2;

        break;

        case 2:

            boundaryGPs1 = gausspoints1;
            boundaryGWs1 = gaussweights1;

            boundaryGPs2 = gpsLeft;
            boundaryGWs2 = gwsLeft;

        break;

        case 3:

            boundaryGPs1 = gausspoints1;
            boundaryGWs1 = gaussweights1;

            boundaryGPs2 = gpsRight;
            boundaryGWs2 = gwsRight;

        break;

        default :
            cout << " Invalid 'side' value in GeomDataHBSplines::getBoundaryGPs2D " << endl;
        break;
    } //switch(side)

  return;
}



void  GeomDataHBSplines::setBoundaryGPs3D(int side, vector<double>& boundaryGPs1, vector<double>& boundaryGWs1, vector<double>& boundaryGPs2, vector<double>& boundaryGWs2, vector<double>& boundaryGPs3, vector<double>& boundaryGWs3)
{
    switch(side)
    {
        case 0:

            boundaryGPs1 = gpsLeft;
            boundaryGWs1 = gwsLeft;

            boundaryGPs2 = gausspoints2;
            boundaryGWs2 = gaussweights2;

            boundaryGPs3 = gausspoints3;
            boundaryGWs3 = gaussweights3;

        break;

        case 1:

            boundaryGPs1 = gpsRight;
            boundaryGWs1 = gwsRight;

            boundaryGPs2 = gausspoints2;
            boundaryGWs2 = gaussweights2;

            boundaryGPs3 = gausspoints3;
            boundaryGWs3 = gaussweights3;

        break;

        case 2:

            boundaryGPs1 = gausspoints1;
            boundaryGWs1 = gaussweights1;

            boundaryGPs2 = gpsLeft;
            boundaryGWs2 = gwsLeft;

            boundaryGPs3 = gausspoints3;
            boundaryGWs3 = gaussweights3;

        break;

        case 3:

            boundaryGPs1 = gausspoints1;
            boundaryGWs1 = gaussweights1;

            boundaryGPs2 = gpsRight;
            boundaryGWs2 = gwsRight;

            boundaryGPs3 = gausspoints3;
            boundaryGWs3 = gaussweights3;

        break;

        case 4:

            boundaryGPs1 = gausspoints1;
            boundaryGWs1 = gaussweights1;

            boundaryGPs2 = gausspoints2;
            boundaryGWs2 = gaussweights2;

            boundaryGPs3 = gpsLeft;
            boundaryGWs3 = gwsLeft;

        break;

        case 5:

            boundaryGPs1 = gausspoints1;
            boundaryGWs1 = gaussweights1;

            boundaryGPs2 = gausspoints2;
            boundaryGWs2 = gaussweights2;

            boundaryGPs3 = gpsRight;
            boundaryGWs3 = gwsRight;

        break;

        default :
            cout << " Invalid 'side' value in GeomDataHBSplines::getBoundaryGPs3D " << endl;
        break;
    } //switch(side)

  return;
}




void  GeomDataHBSplines::setBoundaryGPs1D(int side, GaussQuadrature& quadTemp )
{
    return;
}



void  GeomDataHBSplines::setBoundaryGPs2D(int side, GaussQuadrature& quadTemp)
{
    int  ii=0, jj=0, kk=0;
    myPoint  ptTemp;
    double  wt=0.0;
    ptTemp[2] = 0.0;

    switch(side)
    {
        case 0:

              for(jj=0; jj<gausspoints2.size(); jj++)
              {
                ptTemp[1] = gausspoints2[jj];

                for(ii=0; ii<gpsLeft.size(); ii++)
                {
                  ptTemp[0] = gpsLeft[ii];

                  wt = gaussweights2[jj]*gwsLeft[ii];

                  quadTemp.gausspoints.push_back(ptTemp);
                  quadTemp.gaussweights.push_back(wt);
                }
              }

        break;

        case 1:

              for(jj=0; jj<gausspoints2.size(); jj++)
              {
                ptTemp[1] = gausspoints2[jj];

                for(ii=0; ii<gpsRight.size(); ii++)
                {
                  ptTemp[0] = gpsRight[ii];

                  wt = gaussweights2[jj]*gwsLeft[ii];

                  quadTemp.gausspoints.push_back(ptTemp);
                  quadTemp.gaussweights.push_back(wt);
                }
              }

        break;

        case 2:

              for(jj=0; jj<gpsLeft.size(); jj++)
              {
                ptTemp[1] = gpsLeft[jj];

                for(ii=0; ii<gausspoints1.size(); ii++)
                {
                  ptTemp[0] = gausspoints1[ii];

                  wt = gwsLeft[jj]*gaussweights1[ii];

                  quadTemp.gausspoints.push_back(ptTemp);
                  quadTemp.gaussweights.push_back(wt);
                }
              }

        break;

        case 3:

              for(jj=0; jj<gpsRight.size(); jj++)
              {
                ptTemp[1] = gpsRight[jj];

                for(ii=0; ii<gausspoints1.size(); ii++)
                {
                  ptTemp[0] = gausspoints1[ii];

                  wt = gwsRight[jj]*gaussweights1[ii];

                  quadTemp.gausspoints.push_back(ptTemp);
                  quadTemp.gaussweights.push_back(wt);
                }
              }
        break;

        default :
            cout << " Invalid 'side' value in GeomDataHBSplines::getBoundaryGPs2D " << endl;
        break;
    } //switch(side)

  return;
}



void  GeomDataHBSplines::setBoundaryGPs3D(int side, GaussQuadrature& quadTemp)
{
    int  ii=0, jj=0, kk=0;
    myPoint  ptTemp;
    double  wt=0.0;

    switch(side)
    {
        case 0:

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
                  
                  quadTemp.gausspoints.push_back(ptTemp);
                  quadTemp.gaussweights.push_back(wt);
                }
              }
            }

        break;

        case 1:

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

                  quadTemp.gausspoints.push_back(ptTemp);
                  quadTemp.gaussweights.push_back(wt);
                }
              }
            }

        break;

        case 2:

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

                  quadTemp.gausspoints.push_back(ptTemp);
                  quadTemp.gaussweights.push_back(wt);
                }
              }
            }

        break;

        case 3:

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

                  quadTemp.gausspoints.push_back(ptTemp);
                  quadTemp.gaussweights.push_back(wt);
                }
              }
            }

        break;

        case 4:

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

                  quadTemp.gausspoints.push_back(ptTemp);
                  quadTemp.gaussweights.push_back(wt);
                }
              }
            }

        break;

        case 5:

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

                  quadTemp.gausspoints.push_back(ptTemp);
                  quadTemp.gaussweights.push_back(wt);
                }
              }
            }

        break;

        default :
            cout << " Invalid 'side' value in GeomDataHBSplines::getBoundaryGPs3D " << endl;
        break;
    } //switch(side)

  return;
}




