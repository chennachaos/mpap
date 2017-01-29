#include <Eigen/Dense>

#include <math.h>
#include "Debug.h"
//#include "FunctionsElement.h"
#include "MpapTime.h"
#include "NurbsElem3DStructMixed2fieldStabilised.h"
#include "NurbsShapeFunctions.h"
#include <assert.h>
#include "ComputerTime.h"
#include "TimeFunction.h"
#include "PlotVTK.h"

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
#include <vtkHexahedron.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataSetMapper.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkProperty.h>
#include <vtkFloatArray.h>
#include <vtkLine.h>
#include <vtkLookupTable.h>

#include "util.h"


using namespace std;

extern MpapTime mpapTime;
extern ComputerTime       computerTime;
extern List<TimeFunction> timeFunction;
extern PlotVTK plotvtk;


NurbsElem3DStructMixed2fieldStabilised::NurbsElem3DStructMixed2fieldStabilised(void)
{
  if (debug) cout << " constructor NurbsElem3DStructMixed2fieldStabilised\n\n";

//   cout << " constructor NurbsElem3DStructMixed2fieldStabilised\n\n";

  calcExtraMatrices = true;
}



NurbsElem3DStructMixed2fieldStabilised::~NurbsElem3DStructMixed2fieldStabilised()
{
  if (debug) cout << " destructor NurbsElem3DStructMixed2fieldStabilised\n\n";

  // cout << " destructor NurbsElem3DStructMixed2fieldStabilised\n\n";

}

int NurbsElem3DStructMixed2fieldStabilised::calcStiffnessMatrix(double dt)
{

  return 0;
}



void NurbsElem3DStructMixed2fieldStabilised::contourplot(int index, int nCol, double umin, double umax)
{
  return;
}

int NurbsElem3DStructMixed2fieldStabilised::calcMassMatrix(int lumpInd, double dt)
{
  return 0;
}




int NurbsElem3DStructMixed2fieldStabilised::calcStiffnessAndResidual()
{
  if(finite)
    return  NurbsElem3DStructMixed2fieldStabilised::calcStiffnessAndResidual2();
  else
    return  NurbsElem3DStructMixed2fieldStabilised::calcStiffnessAndResidual1();
}




int NurbsElem3DStructMixed2fieldStabilised::calcStiffnessAndResidual1()
{
   int  err, isw, count, count1, ll, ii, jj, kk, index, gp1, gp2, gp3, mm;

   int  TI, TIp1, TIp2, TIp3, TJ, TJp1, TJp2, TJp3;

   double  F[9], detF, fact, dvol0, dt, Jac, pres, utemp, vtemp, wtemp;
   double  volstr, fact1, fact2, tau, h;
   double  bb1, bb2, bb3, bb4, cc1, cc2, cc3, cc4, dp[3], force[3];

   double  cc[6][6], stre[6], bc[3][6], Idev[6][6], cctmp[6][6];
   double  N[nlbf], dN_dx[nlbf], dN_dy[nlbf],  dN_dz[nlbf];

   Idev3D(Idev);

   double *gaussweights = &(solid0->gaussweights[0]);

  double  *values1 = &(solid1->Values[0][0]);
  double  *values2 = &(solid1->Values[1][0]);
  double  *values3 = &(solid1->Values[2][0]);
  double  *values4 = &(solid1->Values[3][0]);

  int *tt = &(solid0->IEN[elenum][0]);

  double BULK = matDat[0];
  double mu   = matDat[1];
  double eps  = 1.0/BULK;
  //eps = 0.0;
  
  //tau = sqrt(volume());
  
  h = pow(6.0*volume()/PI, 1.0/3.0);

  //h2 = 6.0*volume/PI;
  tau = h*h/(12.0*mu);
  //stabParam /= (degree[0]/degree[0]);

  tau *= elmDat[8];
   //cout << " eps  = " << eps << '\t' << tau << endl;

   Klocal.setZero();
   Flocal.setZero();

   dt = mpapTime.dt;
   count  = 1;
   ll     = 0;
   err    = 0;
   isw    = 3;

   count1 = 0;
   for(gp3=0;gp3<nGP3;gp3++)
   {
   for(gp2=0;gp2<nGP2;gp2++)
   {
   for(gp1=0;gp1<nGP1;gp1++)
   {
          index = count1*3;

          utemp = knotsAtGPs[index];
          vtemp = knotsAtGPs[index+1];
          wtemp = knotsAtGPs[index+2];
          
          solid0->ShapeFunDerivatives(&startindex[0], &(knotsAtGPs[index]), N, dN_dx, dN_dy, dN_dz, Jac);
          
          dvol0 = Jac * gaussweights[count1] * JacMultFact;

          solid1->deformationGradient(&startindex[0], 1, dN_dx, dN_dy, dN_dz, F, detF);

          //F[0] = F[1] = F[2] = 0.0;
          //F[3] = F[4] = F[5] = 0.0;
          //F[6] = F[7] = F[8] = 0.0;
          pres = 0.0;
          dp[0] = dp[1] = dp[2] = 0.0;

          for(ii=0;ii<nlbf;ii++)
          {
            index = tt[ii];

            bb1 = values1[index];
            bb2 = values2[index];
            bb3 = values3[index];
            bb4 = values4[index];

            pres  += N[ii]*bb4;

            dp[0] += dN_dx[ii]*bb4;
            dp[1] += dN_dy[ii]*bb4;
            dp[2] += dN_dz[ii]*bb4;

            //F[0] += bb1*dN_dx[ii];
            //F[1] += bb2*dN_dx[ii];
            //F[2] += bb3*dN_dx[ii];
            //F[3] += bb1*dN_dy[ii];
            //F[4] += bb2*dN_dy[ii];
            //F[5] += bb3*dN_dy[ii];
            //F[6] += bb1*dN_dz[ii];
            //F[7] += bb2*dN_dz[ii];
            //F[8] += bb3*dN_dz[ii];
          }
          //F[0] += 1.0;
          //F[4] += 1.0;
          //F[8] += 1.0;

          //detF = F[0]*(F[4]*F[8] - F[5]*F[7]) - F[3]*(F[1]*F[8] - F[2]*F[7]) + F[6]*(F[1]*F[5] - F[2]*F[4]);

          volstr = (F[0]+F[4]+F[8]-3.0);

          //cout  << volstr << '\t' << Nbar[0] << '\t' << dvol0 << endl;

          matlib3d_(matDat, F, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &isw, &err, &count, (Element*) this);
          if(err !=0)           return 1;

          fact = pres - (stre[0]+stre[1]+stre[2])/3.0 ;

          stre[0] += fact;
          stre[1] += fact;
          stre[2] += fact;

//        printf(" stresses ");        printf("\t%12.8f\t%12.8f\t%12.8f\t%12.8f\n\n", stre[0], stre[1], stre[2], pres);

          for(ii=0;ii<6;ii++)
          {
             for(jj=0;jj<6;jj++)
             {
                cctmp[ii][jj] = 0.0;
                for(mm=0;mm<6;mm++)
                   cctmp[ii][jj] += Idev[ii][mm]*cc[mm][jj];
             }
          }

          for(ii=0;ii<6;ii++)
          {
             for(jj=0;jj<6;jj++)
             {
               cc[ii][jj] = 0.0;
               for(mm=0;mm<6;mm++)
                  cc[ii][jj] += cctmp[ii][mm]*Idev[mm][jj];
             }
          }

          //==============================================
          // CALCULATE TANGENT STIFFNESS and RESIDUAL
          //==============================================

            force[0] = 0.0;
            force[1] = 0.0;
            force[2] = 0.0;

            volstr = volstr - pres*eps;

            for(ii=0;ii<nlbf;ii++)
            {
                bb1 = dN_dx[ii]*dvol0;
                bb2 = dN_dy[ii]*dvol0;
                bb3 = dN_dz[ii]*dvol0;
                bb4 = N[ii]*dvol0;
                
                for(kk=0;kk<6;kk++)
                {
                  bc[0][kk] = bb1 * cc[0][kk] + bb2 * cc[3][kk] + bb3 * cc[5][kk];
                  bc[1][kk] = bb1 * cc[3][kk] + bb2 * cc[1][kk] + bb3 * cc[4][kk];
                  bc[2][kk] = bb1 * cc[5][kk] + bb2 * cc[4][kk] + bb3 * cc[2][kk];
                }

                TI   = 4*ii;
                TIp1 = TI+1;
                TIp2 = TI+2;
                TIp3 = TI+3;

                Flocal[TI]   += (bb4*force[0] - bb1*stre[0] - bb2*stre[3] - bb3*stre[5]) ;
                Flocal[TIp1] += (bb4*force[1] - bb1*stre[3] - bb2*stre[1] - bb3*stre[4]) ;
                Flocal[TIp2] += (bb4*force[2] - bb1*stre[5] - bb2*stre[4] - bb3*stre[2]) ;
                Flocal[TIp3] -= (bb4*volstr);

                // PSPG stabilization
                Flocal[TIp3] += tau*(bb1*dp[0] + bb2*dp[1] + bb3*dp[2]);

                for(jj=0;jj<nlbf;jj++)
                {
                   cc1 = dN_dx[jj];
                   cc2 = dN_dy[jj];
                   cc3 = dN_dz[jj];
                   cc4 = N[jj];

                   TJ   = 4*jj;
                   TJp1 = TJ+1;
                   TJp2 = TJ+2;
                   TJp3 = TJ+3;

                   Klocal(TI,   TJ)    +=  (bc[0][0] * cc1 + bc[0][3] * cc2 + bc[0][5] * cc3) ;
                   Klocal(TI,   TJp1)  +=  (bc[0][1] * cc2 + bc[0][3] * cc1 + bc[0][4] * cc3) ;
                   Klocal(TI,   TJp2)  +=  (bc[0][2] * cc3 + bc[0][4] * cc2 + bc[0][5] * cc1) ;
                   Klocal(TI,   TJp3)  +=  (bb1 * cc4);

                   Klocal(TIp1, TJ)    +=  (bc[1][0] * cc1 + bc[1][3] * cc2 + bc[1][5] * cc3) ;
                   Klocal(TIp1, TJp1)  +=  (bc[1][1] * cc2 + bc[1][3] * cc1 + bc[1][4] * cc3) ;
                   Klocal(TIp1, TJp2)  +=  (bc[1][2] * cc3 + bc[1][4] * cc2 + bc[1][5] * cc1) ;
                   Klocal(TIp1, TJp3)  +=  (bb2 * cc4);

                   Klocal(TIp2, TJ)    +=  (bc[2][0] * cc1 + bc[2][3] * cc2 + bc[2][5] * cc3) ;
                   Klocal(TIp2, TJp1)  +=  (bc[2][1] * cc2 + bc[2][3] * cc1 + bc[2][4] * cc3) ;
                   Klocal(TIp2, TJp2)  +=  (bc[2][2] * cc3 + bc[2][4] * cc2 + bc[2][5] * cc1) ;
                   Klocal(TIp2, TJp3)  +=  (bb3 * cc4);

                   Klocal(TIp3, TJ)    +=  (bb4 * cc1);
                   Klocal(TIp3, TJp1)  +=  (bb4 * cc2);
                   Klocal(TIp3, TJp2)  +=  (bb4 * cc3);
                   Klocal(TIp3, TJp3)  -=  (bb4 * cc4)*eps;

                   // PSPG stabilization
                   Klocal(TIp3, TJp3)   -= tau*(bb1*cc1 + bb2*cc2 + bb3*cc3);
                }
            }

          count++;
          count1++;
          ll += nivGP;
   } // gp1
   } // gp2
   } // gp3

  return 0;
}





int NurbsElem3DStructMixed2fieldStabilised::calcStiffnessAndResidual2()
{
   int  err, isw, count, count1, ll, ii, jj, kk, index, gp1, gp2, gp3, mm;
   int  TI, TIp1, TIp2, TIp3, TJ, TJp1, TJp2, TJp3;

   double  F[9], detF, dvol0, dt, Jac, dummy, pres, utemp, vtemp, wtemp, r1d3 = 1.0/3.0, r2d3 = 2.0*r1d3;
   double  fact, fact1, fact2, fact3, pbar, dvol, tau, h;
   double  bb1, bb2, bb3, bb4, bb5, cc1, cc2, cc3, cc4, volstr;

   double  D11[6][6], stre[6], bc[3][6], Idev[6][6], cctmp[6][6], strdev[6], force[3], dp[3];
   double  N[nlbf], dN_dx[nlbf], dN_dy[nlbf],  dN_dz[nlbf];

   Idev3D(Idev);

   double *gaussweights = &(solid0->gaussweights[0]);

  double  *values1 = &(solid1->Values[0][0]);
  double  *values2 = &(solid1->Values[1][0]);
  double  *values3 = &(solid1->Values[2][0]);
  double  *values4 = &(solid1->Values[3][0]);

  int *tt = &(solid0->IEN[elenum][0]);


   double BULK = matDat[0];
   double mu   = matDat[1];
   double eps  = 1.0/BULK;
   //eps = 0.0;

  h = pow(3.0*volume()/PI/4.0, 1.0/3.0);

  //h = 6.0*volume()/PI;
  tau = h*h/(12.0*mu);
  //stabParam /= (degree[0]/degree[0]);

  tau *= elmDat[7];
   //cout << " eps  = " << eps << '\t' << tau << endl;

   Klocal.setZero();
   Flocal.setZero();

   count  = 1;
   ll     = 0;
   err    = 0;
   isw    = 3;
   dt = mpapTime.dt;

   count1 = 0;
   for(gp3=0;gp3<nGP3;gp3++)
   {
   for(gp2=0;gp2<nGP2;gp2++)
   {
   for(gp1=0;gp1<nGP1;gp1++)
   {
          index = 3*count1;

          utemp = knotsAtGPs[index];
          vtemp = knotsAtGPs[index+1];
          wtemp = knotsAtGPs[index+2];

          solid0->ShapeFunDerivatives(&startindex[0], &(knotsAtGPs[index]), N, dN_dx, dN_dy, dN_dz, Jac);

          dvol0 = Jac * gaussweights[count1] * JacMultFact;

          solid1->deformationGradient(&startindex[0], 1, dN_dx, dN_dy, dN_dz, F, detF);

          index = 3*count1;
          solid1->ShapeFunDerivatives(&startindex[0], &(knotsAtGPs[index]), N, dN_dx, dN_dy, dN_dz, Jac);

          dvol = Jac * gaussweights[count1] * JacMultFact;

          //F[0] = F[1] = F[2] = 0.0;
          //F[3] = F[4] = F[5] = 0.0;
          //F[6] = F[7] = F[8] = 0.0;
          pres = 0.0;
          dp[0] = dp[1] = dp[2] = 0.0;

          for(ii=0;ii<nlbf;ii++)
          {
            index = tt[ii];

            bb1 = values1[index];
            bb2 = values2[index];
            bb3 = values3[index];
            bb4 = values4[index];

            pres  += N[ii]*bb4;

            dp[0] += dN_dx[ii]*bb4;
            dp[1] += dN_dy[ii]*bb4;
            dp[2] += dN_dz[ii]*bb4;

            //F[0] += bb1*dN_dx[ii];
            //F[1] += bb2*dN_dx[ii];
            //F[2] += bb3*dN_dx[ii];
            //F[3] += bb1*dN_dy[ii];
            //F[4] += bb2*dN_dy[ii];
            //F[5] += bb3*dN_dy[ii];
            //F[6] += bb1*dN_dz[ii];
            //F[7] += bb2*dN_dz[ii];
            //F[8] += bb3*dN_dz[ii];
          }
          //F[0] += 1.0;
          //F[4] += 1.0;
          //F[8] += 1.0;

          //detF = F[0]*(F[4]*F[8] - F[5]*F[7]) - F[3]*(F[1]*F[8] - F[2]*F[7]) + F[6]*(F[1]*F[5] - F[2]*F[4]);

          //volstr = (F[0]+F[4]+F[8]-3.0);


          matlib3d_(matDat, F, stre, D11[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &isw, &err, &count, (Element*) this);
          if(err !=0)           return 1;

          //cout << " pres " << pres << endl;

          pbar = (stre[0]+stre[1]+stre[2])/3.0;

          strdev[0] = stre[0] - pbar;
          strdev[1] = stre[1] - pbar;
          strdev[2] = stre[2] - pbar;
          strdev[3] = stre[3];
          strdev[4] = stre[4];
          strdev[5] = stre[5];

          stre[0] = strdev[0] + pres;
          stre[1] = strdev[1] + pres;
          stre[2] = strdev[2] + pres;

//        printf(" stresses ");        printf("\t%12.8f\t%12.8f\t%12.8f\t%12.8f\n\n", stre[0], stre[1], stre[2], pres);
//        printf(" values   ");        printf("\t%12.8f\t%12.8f\t%12.8f\t%12.8f\n\n", detF, detFbar, pbar, dvol);

          for(ii=0;ii<6;ii++)
          {
             for(jj=0;jj<6;jj++)
             {
                cctmp[ii][jj] = 0.0;
                for(mm=0;mm<6;mm++)
                  cctmp[ii][jj] += Idev[ii][mm] * D11[mm][jj];
             }
          }

          for(ii=0;ii<6;ii++)
          {
             for(jj=0;jj<6;jj++)
             {
                D11[ii][jj] = 0.0;
                for(mm=0;mm<6;mm++)
                  D11[ii][jj] += cctmp[ii][mm] * Idev[mm][jj];
             }
          }

          fact = (pbar - pres);
          fact1 = (r2d3*pbar - pres);
          fact2 = 2.0 * fact;

          for(ii=0;ii<3;ii++)
          {
             for(jj=0;jj<6;jj++)
             {
                D11[ii][jj] -= r2d3*strdev[jj] ;
                D11[jj][ii] -= r2d3*strdev[jj] ;
             }
             for(jj=0;jj<3;jj++)
                D11[ii][jj] -= fact1;

             D11[ii][ii] += fact2;
             index = ii+3;
             D11[index][index] += fact;
          }

          //==============================================
          // CALCULATE TANGENT STIFFNESS and RESIDUAL
          //==============================================
            
            force[0] = 0.0;//1*timeFunction[0].prop;
            force[1] = 0.0;
            force[2] = 0.0;

            volstr  = detF - 1.0 - pres*eps;

            for(ii=0;ii<nlbf;ii++)
            {
                bb1 = dN_dx[ii]*dvol;
                bb2 = dN_dy[ii]*dvol;
                bb3 = dN_dz[ii]*dvol;
                bb4 = N[ii]*dvol;
                bb5 = N[ii]*dvol0;

                for(kk=0;kk<6;kk++)
                {
                  bc[0][kk] = (bb1 * D11[0][kk] + bb2 * D11[3][kk] + bb3 * D11[5][kk]);
                  bc[1][kk] = (bb1 * D11[3][kk] + bb2 * D11[1][kk] + bb3 * D11[4][kk]);
                  bc[2][kk] = (bb1 * D11[5][kk] + bb2 * D11[4][kk] + bb3 * D11[2][kk]);
                }

                TI   = 4*ii;
                TIp1 = TI+1;
                TIp2 = TI+2;
                TIp3 = TI+3;

                fact1 = bb1*stre[0] + bb2*stre[3] + bb3*stre[5] ;
                fact2 = bb1*stre[3] + bb2*stre[1] + bb3*stre[4] ;
                fact3 = bb1*stre[5] + bb2*stre[4] + bb3*stre[2] ;

                Flocal[TI]   += (bb5*force[0] - fact1) ;
                Flocal[TIp1] += (bb5*force[1] - fact2) ;
                Flocal[TIp2] += (bb5*force[2] - fact3) ;
                Flocal[TIp3] -=  bb5*volstr ;

                // PSPG stabilization
                Flocal[TIp3] += tau*(bb1*dp[0] + bb2*dp[1] + bb3*dp[2]);


                for(jj=0;jj<nlbf;jj++)
                {
                   cc1 = dN_dx[jj];
                   cc2 = dN_dy[jj];
                   cc3 = dN_dz[jj];
                   cc4 = N[jj];

                   TJ   = 4*jj;
                   TJp1 = TJ+1;
                   TJp2 = TJ+2;
                   TJp3 = TJ+3;

                   fact = fact1 * bb1 + fact2 * bb2 + fact3 * bb3;

                   Klocal(TI,   TJ)    +=  (bc[0][0] * cc1 + bc[0][3] * cc2 + bc[0][5] * cc3 + fact) ;
                   Klocal(TI,   TJp1)  +=  (bc[0][1] * cc2 + bc[0][3] * cc1 + bc[0][4] * cc3) ;
                   Klocal(TI,   TJp2)  +=  (bc[0][2] * cc3 + bc[0][4] * cc2 + bc[0][5] * cc1) ;
                   Klocal(TI,   TJp3)  +=  (bb1 * cc4);

                   Klocal(TIp1, TJ)    +=  (bc[1][0] * cc1 + bc[1][3] * cc2 + bc[1][5] * cc3) ;
                   Klocal(TIp1, TJp1)  +=  (bc[1][1] * cc2 + bc[1][3] * cc1 + bc[1][4] * cc3 + fact) ;
                   Klocal(TIp1, TJp2)  +=  (bc[1][2] * cc3 + bc[1][4] * cc2 + bc[1][5] * cc1) ;
                   Klocal(TIp1, TJp3)  +=  (bb2 * cc4);

                   Klocal(TIp2, TJ)    +=  (bc[2][0] * cc1 + bc[2][3] * cc2 + bc[2][5] * cc3) ;
                   Klocal(TIp2, TJp1)  +=  (bc[2][1] * cc2 + bc[2][3] * cc1 + bc[2][4] * cc3) ;
                   Klocal(TIp2, TJp2)  +=  (bc[2][2] * cc3 + bc[2][4] * cc2 + bc[2][5] * cc1 + fact) ;
                   Klocal(TIp2, TJp3)  +=  (bb3 * cc4);

                   Klocal(TIp3, TJ)    +=  (bb4 * cc1);
                   Klocal(TIp3, TJp1)  +=  (bb4 * cc2);
                   Klocal(TIp3, TJp2)  +=  (bb4 * cc3);
                   Klocal(TIp3, TJp3)  -=  (bb5 * cc4)*eps;

                   // PSPG stabilization
                   Klocal(TIp3, TJp3)   -= tau*(bb1*cc1 + bb2*cc2 + bb3*cc3);
                }
            }

          count++;
          count1++;
          ll += nivGP;
   } // gp1
   } // gp2
   } // gp3

  //printMatrix(Klocal); printf("\n\n"); printVector(Flocal);

  return 0;
}





int NurbsElem3DStructMixed2fieldStabilised::calcInternalForces()
{
/*
   double BULK = matDat[0];

   int  err, isw, count, count1, ll, ii, jj, kk, index, gp1, gp2, gp3, mm;

   int  threeI, threeIp1, threeIp2;

   int  sizep = solid2->nlbf;

   double  F[9], detF, dvol0, dt, Jac, dummy, pres, bb1, bb2, bb3, utemp, vtemp, wtemp, r1d3 = 1.0/3.0, r2d3 = 2.0*r1d3;
   
   double  fact, fact1, fact2, fact3, pbar, dvol;

   double  D11[6][6], stre[6], N[nlbf], dN_dx[nlbf], dN_dy[nlbf],  dN_dz[nlbf], Nbar[sizep];

   double *gaussweights = &(solid0->gaussweights[0]);

   resi.zero();
   resi2.zero();

   count  = 1;
   ll     = 0;
   err    = 0;
   isw    = 3;
   dt = mpapTime.dt;

   count1 = 0;
   for(gp3=0;gp3<nGP3;gp3++)
   {
   for(gp2=0;gp2<nGP2;gp2++)
   {
   for(gp1=0;gp1<nGP1;gp1++)
   {
          index = 3*count1;

          utemp = knotsAtGPs[index];
          vtemp = knotsAtGPs[index+1];
          wtemp = knotsAtGPs[index+2];

          solid1->ShapeFunDerivatives(&startindex[0], &(knotsAtGPs[index]), N, dN_dx, dN_dy, dN_dz, Jac);

          dvol = Jac * gaussweights[count1] * JacMultFact;

          solid0->deformationGradient(&startindex[0], 0, dN_dx, dN_dy, dN_dz, F, detF);

          dvol0 = dvol/detF;

          matlib3d_(matDat, F, stre, D11[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &isw, &err, &count, (Element*) this);
          if(err !=0)           return 1;

          pres = solid2->computeValueAndShanpeFns(1, utemp, vtemp, wtemp, Nbar);

          pbar = (stre[0]+stre[1]+stre[2])/3.0;
          fact = pbar - pres;

          stre[0] -= fact;
          stre[1] -= fact;
          stre[2] -= fact;

//        printf(" stresses ");        printf("\t%12.8f\t%12.8f\t%12.8f\t%12.8f\n\n", stre[0], stre[1], stre[2], pres);
//        printf(" values   ");        printf("\t%12.8f\t%12.8f\t%12.8f\t%12.8f\n\n", detF, detFbar, pbar, dvol);

          for(ii=0;ii<6;ii++)
             stre[ii] *= dvol;

          for(ii=0;ii<nlbf;ii++)
          {
                bb1 = dN_dx[ii];
                bb2 = dN_dy[ii];
                bb3 = dN_dz[ii];
                
                threeI   = 3*ii;
                threeIp1 = threeI+1;
                threeIp2 = threeI+2;

                resi[threeI]   -= (bb1*stre[0] + bb2*stre[3] + bb3*stre[5]);
                resi[threeIp1] -= (bb1*stre[3] + bb2*stre[1] + bb3*stre[4]) ;
                resi[threeIp2] -= (bb1*stre[5] + bb2*stre[4] + bb3*stre[2]) ;
          }

            fact1 = (detF - 1.0 - pres/BULK) * dvol0;
          
            for(ii=0;ii<sizep;ii++)
                resi2[ii] -= Nbar[ii] * fact1;

          count++;
          count1++;
          ll += nivGP;
   } // gp1
   } // gp2
   } // gp3
*/
//  printForceVector();

  return 0;
}




int NurbsElem3DStructMixed2fieldStabilised::calcOutput(double u1, double v1)
{
  return 0;
}




void NurbsElem3DStructMixed2fieldStabilised::discreteContourplot(int vartype, int varindex, int index, int nCol, double umin, double umax)
{
   if(index > nivGP)
   {
      cout << '\t' << " Error in NurbsElem3DStructSolid::contourplot " << endl;
      return;
   }
  
   double outval[500];

   switch(vartype)
   {
       case 0:  // plot total strain
       case 1:  // plot elastic strain
       case 2:  // plot plastic strain

                projectStrain(vartype, varindex, outval);

              break;

       case 3:  // plot stress

                projectStress(varindex, outval);

              break;

       case 4:  // plot element internal variables


                projectIntVar(index, outval);

              break;

       default:

              cout  << "           Invalid Variable Type to project " << endl;
              break;

   }

    double  uu, vv, ww, du, dv, dw, utemp, vtemp, wtemp;
    int count=0, ii, kk, jj, ll, n1, n2, n3, nn, ind1, ind2, ind3, ind4, ind5, ind6;

/*
  count=0;
  for(kk=0;kk<nGP3;kk++)
  {
     for(jj=0;jj<nGP2;jj++)
     { 
        for(ii=0;ii<nGP1;ii++)
        {
           printf("\t%12.6f",outval[count++]);
        }
        printf("\n");
     }
     printf("\n");
     printf("\n");
  }

  cout << endl;
  cout << endl;
*/
    du = uvalues[2]/nGP1;
    dv = vvalues[2]/nGP2;
    dw = wvalues[2]/nGP3;

    n1 = nGP1+1;
    n2 = nGP2+1;
    n3 = nGP3+1;

    vtkSmartPointer<vtkPoints>            points    =  vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkUnstructuredGrid>  uGrid     =  vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkHexahedron>        hex       =  vtkSmartPointer<vtkHexahedron>::New();
    vtkSmartPointer<vtkDataSetMapper>     mapper1   =  vtkSmartPointer<vtkDataSetMapper>::New();
    vtkSmartPointer<vtkActor>             actor1    =  vtkSmartPointer<vtkActor>::New();
    vtkSmartPointer<vtkLookupTable>       hueLut    =  vtkSmartPointer<vtkLookupTable>::New();
    vtkSmartPointer<vtkFloatArray>        scalars   =  vtkSmartPointer<vtkFloatArray>::New();


    vtkIdType pts[8], pt1, pt2, start;
    
    start = plotvtk.points->GetNumberOfPoints();

    EPOINT  EP1;

    nn = (nGP1+1)*(nGP2+1);

    scalars->SetNumberOfTuples(nn*n3);


    // Calculate the Points on the solid
    count = 0;
    if(finite)
    {
        for(kk=0;kk<n3;kk++)
        {
            wtemp = wvalues[0]+dw*kk;
            for(jj=0;jj<n2;jj++)
            {
               vtemp = vvalues[0]+dv*jj;
               for(ii=0;ii<n1;ii++)
               {
                  EP1 = solid1->SolidPoint(uvalues[0]+ii*du, vtemp, wtemp).CalcEuclid();
                  points->InsertNextPoint(EP1.x, EP1.y, EP1.z);

                  //scalars->SetTuple1(count, outval[count]);
                  //count++;
               }
            }
        }
    }
    else
    {
        for(kk=0;kk<n3;kk++)
        {
            wtemp = wvalues[0]+dw*kk;
            for(jj=0;jj<n2;jj++)
            {
               vtemp = vvalues[0]+dv*jj;
               for(ii=0;ii<n1;ii++)
               {
                  EP1 = solid0->SolidPoint(uvalues[0]+ii*du, vtemp, wtemp).CalcEuclid();
                  points->InsertNextPoint(EP1.x, EP1.y, EP1.z);

                  //scalars->SetTuple1(count, outval[count]);
                  //count++;
               }
            }
        }
    }

    uGrid->SetPoints(points);

    count=0;
    for(kk=0;kk<nGP3;kk++)
    {
        ind5 = nn*kk;
        ind6 = nn*(kk+1);

        for(jj=0;jj<nGP2;jj++)
        {
           ind1 = ind5 + n1*jj;
           ind2 = ind5 + n1*(jj+1);
       
           ind3 = ind6 + n1*jj;
           ind4 = ind6 + n1*(jj+1);

           for(ii=0;ii<nGP1;ii++)
           {
              pts[0] = ind1+ii;          pts[4] = ind3+ii;

              pts[1] = pts[0]+1;         pts[5] = pts[4]+1;
              pts[3] = ind2+ii;          pts[7] = ind4+ii;
              pts[2] = pts[3]+1;         pts[6] = pts[7]+1;
              
              for(ll=0;ll<8;ll++)
              {
                 scalars->SetTuple1(pts[ll], outval[count]);
                 hex->GetPointIds()->SetId(ll, pts[ll]);
              }
          
              uGrid->InsertNextCell(hex->GetCellType(), hex->GetPointIds());
              count++;
           }
        }
    }

    uGrid->GetPointData()->SetScalars(scalars);

#if VTK_MAJOR_VERSION == 5
    mapper1->SetInputConnection(uGrid->GetProducerPort());
#else
    mapper1->SetInputData(uGrid);
#endif

    actor1->SetMapper(mapper1);
    
    mapper1->SetScalarRange(umin, umax);

    hueLut->SetHueRange(0.66667, 0.0);
    //hueLut->SetTableRange (0, 3);
    //hueLut->SetNumberOfTableValues(nCol);
    hueLut->SetRampToLinear();
    //hueLut->SetScaleToLinear();
    hueLut->Build();
 
    mapper1->SetLookupTable(hueLut);
    scalars->SetLookupTable(hueLut);

//    actor1->GetProperty()->SetLineWidth(2.0);
//    actor1->GetProperty()->EdgeVisibilityOn();

    plotvtk.rendr->AddActor(actor1);

  return;
}






void NurbsElem3DStructMixed2fieldStabilised::projectToKnots(bool extrapolateFlag, int vartype, int varindex, int index)
{
/*
   vals2project[0] = intVar2[indx];
   vals2project[1] = intVar2[(nGP1-1)*nivGP+indx];
   vals2project[2] = intVar2[nGP1*(nGP2-1)*nivGP+indx];
   vals2project[3] = intVar2[(nGP1*nGP2-1)*nivGP+indx];
*/

   double outval[500];

   switch(vartype)
   {
       case 1:  // plot total strain
       case 2:  // plot elastic strain
       case 3:  // plot plastic strain

                projectStrain(vartype, varindex, outval);

              break;

       case 4:  // plot stress

                projectStress(varindex, outval);

              break;

       case 5:  // plot element internal variables

                projectIntVar(index, outval);

              break;

       default:

              cout  << "           Invalid Variable Type to project " << endl;
              break;

   }


   assert(vals2project.n == 8);

    if(extrapolateFlag)
    {
       for(int ii=0;ii<8;ii++)
         vals2project[ii] = extrapolate3D(nGP1, nGP2, nGP3, (ii+1), outval);
    }
    else
    {
       vals2project[0] = outval[0];
       vals2project[1] = outval[nGP1-1];
       vals2project[2] = outval[nGP1*(nGP2-1)];
       vals2project[3] = outval[nGP1*nGP2-1];
       
       int  nn = nGP1*nGP2*(nGP3-1);

       vals2project[4] = outval[nn];
       vals2project[5] = outval[nn+nGP1-1];
       vals2project[6] = outval[nn+nGP1*(nGP2-1)];
       vals2project[7] = outval[nGP-1];
    }

//   cout << '\t' << vals2project << endl; cout << endl;

  return;
}




void NurbsElem3DStructMixed2fieldStabilised::projectStress(int varindex, double* outval)
{

    if(varindex > 8)
    {
       cout << '\t' << "    NurbsElem3DStructSolid::projectStress .... : Error in 'varindex' " << endl;
       return;
    }

   double F[9], detF, stre[6], cc[6][6], Jac, dt, pres, utemp, pbar;

   int   err    = 0,
         isw    = 3,
         count  = 1, 
         count1 = 0, index, ll = 0, gp1, gp2, gp3;

   double  N[nlbf], dN_dx[nlbf], dN_dy[nlbf], dN_dz[nlbf];

   dt = mpapTime.dt;

   int nivEL = nGP * nivGP;
   for(ll=0;ll<nivEL;ll++)
     intVar2[ll] = intVar1[ll];

  
    ll = 0;
    for(gp3=0;gp3<nGP3;gp3++)
    {
    for(gp2=0;gp2<nGP2;gp2++)
    {
    for(gp1=0;gp1<nGP1;gp1++)
    {
          index = count1*3;

          solid0->ShapeFunDerivatives(&startindex[0], &(knotsAtGPs[index]), N, dN_dx, dN_dy, dN_dz, Jac);

          solid1->deformationGradient(&startindex[0], 1, dN_dx, dN_dy, dN_dz, F, detF);

          matlib3d_(matDat, F, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &isw, &err, &count, NULL);

          pres = solid1->computeValue(4, knotsAtGPs[index], knotsAtGPs[index+1], knotsAtGPs[index+2]);

          pbar = (stre[0]+stre[1]+stre[2])/3.0;

          utemp = pres - pbar;

          stre[0] += utemp ;
          stre[1] += utemp ;
          stre[2] += utemp ;

          if(varindex < 6)
             outval[count1] = stre[varindex];
          else if(varindex == 6)
             outval[count1] = vonMises3D(stre);
          else if(varindex == 7)
             outval[count1] = pres;

          count++;
          count1++;
          ll += nivGP;

    }//gp1
    }//gp2
    }//gp3

  return;
}




void NurbsElem3DStructMixed2fieldStabilised::projectStrain(int vartype, int varindex, double* outval)
{
  return;
}



void NurbsElem3DStructMixed2fieldStabilised::projectIntVar(int index, double* outval)
{
   int ind1, ind2, ii, jj, kk;

   ind1 = 0;
   for(kk=0;kk<nGP3;kk++)
   {
      for(jj=0;jj<nGP2;jj++)
      {
          for(ii=0;ii<nGP1;ii++)
          {
              outval[ind1] = intVar2[ind1*nivGP+index];
              ind1++;
          }
      }
   }

   return;
}




void NurbsElem3DStructMixed2fieldStabilised::AssembleElementMatrix2(int index, MatrixXd& stiff, MatrixXd& CC)
{
    int row, aa, bb, *colC,  ind;

    ind  = surf2->IEN[elenum].n;
    colC = &(surf2->IEN[elenum][0]);

    for(aa=0;aa<nsize;aa++)
    {
        row = forassy[aa];
        for(bb=0;bb<nsize;bb++)
            stiff(row, forassy[bb]) += stiffness_local[aa][bb];

        for(bb=0;bb<ind;bb++)
           CC(row, colC[bb]) += Kup(aa,bb);
    }

  return;
}



void NurbsElem3DStructMixed2fieldStabilised::AssembleElementVector2(bool firstIter, int flag, VectorXd& rhs, double* reac, VectorXd& uu)
{
   // flag == true  ---> just external force vector
   // flag == false ---> internal load vector

   if(flag == 1)
   {
      for(int aa=0;aa<nsize;aa++)
      {
         rhs(forassy[aa]) += Flocal[aa];
      }
   }
   else if(flag == 0)
   {
      for(int aa=0;aa<nsize;aa++)
      {
         rhs(forassy[aa])  += Flocal[aa];
         reac[forassy[aa]] += Flocal[aa];
      }
   }
   else if(flag == 2)
   {
      int  size2 = surf2->nsize;
      int  *tt2  = &(surf2->LM[elenum][0]);
      int aa, bb, ind;

      for(aa=0;aa<size2;aa++)
      {
          ind = tt2[aa];

          for(bb=0;bb<nsize;bb++)
          {
             rhs[ind] += Kup(bb,aa) * uu(forassy[bb]);
          }
      }
   }
   else
   {
      cerr << " ERROR in flag in 'NurbsElem3DStructMixed2fieldStabilised::AssembleElementVector2' " << endl;
      return;
   }


//  cout << " resi " << resi << endl;

  return;
}






void NurbsElem3DStructMixed2fieldStabilised::AssembleElementMatrix(int index, MatrixSparseArray<double>& mtx)
{
     int nn=0, aa, bb, size2;

     size2 = solid2->nsize;

     for(aa=0;aa<nsize;aa++)
     {
         for(bb=0;bb<nsize;bb++)
         {
             nn = forassembly[aa][bb];
             if(nn != -1)
                mtx.x[nn-1] += Klocal(aa, bb);
         }

         for(bb=0;bb<size2;bb++)
         {
             nn = forassyKup[aa][bb];
             if(nn != -1)
                mtx.x[nn-1] += Kup(aa,bb);

             nn = forassyKpu[bb][aa];
             if(nn != -1)
                mtx.x[nn-1] += Kup(aa,bb);
         }
     }
    for(aa=0;aa<size2;aa++)
    {
      for(bb=0;bb<size2;bb++)
      {
        nn = forassyKtt[aa][bb];
        if(nn != -1)
          mtx.x[nn-1] += Kpp(aa,bb);
      }
    }

  return;
}





void NurbsElem3DStructMixed2fieldStabilised::AssembleElementVector(bool firstIter, bool flag, double* rhs, double* reac, int start1, int start2)
{
   // flag == true  ---> just external force vector
   // flag == false ---> internal load vector + contributions from nodes with specified displacement BCs

   int *tt;

   tt = &(solid0->LM[elenum][0]);

   if(flag)
   {
      for(int aa=0;aa<nsize;aa++)
      {
         if(tt[aa] != -1)
            rhs[tt[aa]] += Flocal[aa];
      }
   }
   else
   {
      double fact;
      int aa, bb, ind;

      for(aa=0;aa<nsize;aa++)
      {
         if(tt[aa] != -1)
           rhs[tt[aa]] += Flocal[aa];

         // add up reaction forces
         reac[forassy[aa]] += Flocal[aa];
      }

      // contribution to the pressure variables from the applied displacements

      if(firstIter)
      {
         for(aa=0;aa<nsize;aa++)
         {
             if(tt[aa] == -1)
             {
                 fact = mpapTime.dt * primvar[aa];
                 fact = primvar[aa];
                 for(bb=0;bb<nsize;bb++)
                 {
                    if(tt[bb] != -1)
                      rhs[tt[bb]] -= Klocal(bb,aa) * fact;
                 }
             }
         }
      }
   }

//  cout << " resi " << resi << endl;

  tt= NULL;

  return;
}



/*
void  NurbsElem3DStructMixed2fieldStabilised::AssembleElementMatrix(int index, Mat mtx, int start1, int start2)
{
    PetscErrorCode ierr;
    int  ii, jj, nn=0, aa, bb, size2, ind, *tt1, *tt2;

    tt1 = &(solid0->LM[elenum][0]);
    tt2 = &(solid2->LM[elenum2][0]);
    size2 = solid2->nsize;

    //cout << elenum << '\t' << elenum2 << endl;
    //cout << nsize << '\t' << size2 << endl;
    //cout << surf0->LM[elenum] << endl;
    //cout << surf2->LM[elenum2] << endl;
    
    //printMatrix(Klocal);
    //printf("\n\n");
    //printMatrix(Kup);
    //printf("\n\n");
    //printMatrix(Kpp);
    //printf("\n\n");

    for(ii=0;ii<nsize;ii++)
    {
      aa = tt1[ii];
      if( aa != -1)
      {
        for(jj=0;jj<nsize;jj++)
        {
          bb = tt1[jj];
          if(bb != -1)
          {
            //cout << ii << '\t' << jj << '\t' << aa << '\t' << bb << endl;
            ierr = MatSetValues(mtx, 1, &aa, 1, &bb, &(Klocal(ii, jj)), ADD_VALUES);
          }
        }
        //cout << " ii = " << ii << endl;

        for(jj=0;jj<size2;jj++)
        {
          bb = tt2[jj];
          if(bb != -1)
          {
            bb += start1;
            ierr = MatSetValues(mtx, 1, &aa, 1, &bb, &(Kup(ii, jj)), ADD_VALUES);
            ierr = MatSetValues(mtx, 1, &bb, 1, &aa, &(Kup(ii, jj)), ADD_VALUES);
          }
        }
        //cout << " ii = " << ii << endl;
      }
    }
    //cout << " qqqqqqqqqqqq " << endl;
    for(ii=0;ii<size2;ii++)
    {
      aa = tt2[ii];
      if( aa != -1)
      {
        aa += start1;
        for(jj=0;jj<size2;jj++)
        {
          bb = tt2[jj];
          if(bb != -1)
          {
            bb += start1;
            ierr = MatSetValues(mtx, 1, &aa, 1, &bb, &(Kpp(ii, jj)), ADD_VALUES);
          }
        }
      }
    }

  return;
}
*/



void NurbsElem3DStructMixed2fieldStabilised::toPostprocess(int vartype, int varindex, int type, SparseMatrixXd&  coeffMat, VectorXd& rhsVec)
{
/*
   MatrixXd  Nlocal(nlbf,nlbf);
   VectorXd  NN(nlbf), rhslocal(nlbf);

   Nlocal.setZero();   
   rhslocal.setZero();


   double F[9], detF, stre[6], cc[6][6], Jac, dt, pres, utemp, vtemp, wtemp, pbar;

   int   err    = 0,
         isw    = 3,
         count  = 1, 
         count1 = 0, index, ll = 0, gp1, gp2, gp3, ii, jj, row, col;

   double  N[nlbf], dN_dx[nlbf], dN_dy[nlbf], dN_dz[nlbf];

   dt = mpapTime.dt;

   int nivEL = nGP * nivGP;
   for(ll=0;ll<nivEL;ll++)
     intVar2[ll] = intVar1[ll];

  
    ll = 0;
    for(gp3=0;gp3<nGP3;gp3++)
    {
    for(gp2=0;gp2<nGP2;gp2++)
    {
    for(gp1=0;gp1<nGP1;gp1++)
    {
          index = count1*3;

          utemp = knotsAtGPs[index];
          vtemp = knotsAtGPs[index+1];
          wtemp = knotsAtGPs[index+2];


          solid0->ShapeFunDerivatives(&startindex[0], &(knotsAtGPs[index]), N, dN_dx, dN_dy, dN_dz, Jac);

          solid0->ShapeFunctions(utemp, vtemp, wtemp, &NN(0));

          solid1->deformationGradient(&startindex[0], 1, dN_dx, dN_dy, dN_dz, F, detF);

          matlib3d_(matDat, F, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &isw, &err, &count, NULL);

          pres = solid2->computeValue(1, utemp, vtemp, wtemp);

          pbar = (stre[0]+stre[1]+stre[2])/3.0;
          
          utemp = pres - pbar;

          stre[0] += utemp ;
          stre[1] += utemp ;
          stre[2] += utemp ;

          Nlocal += NN * NN.transpose();
        
          if(varindex < 6)
             rhslocal += NN * stre[varindex];
          else if(varindex == 6)
             rhslocal += NN * vonMises3D(stre);
          else if(varindex == 7)
             rhslocal += NN * pres;

          count++;
          count1++;
          ll += nivGP;

    }//gp1
    }//gp2
    }//gp3


    int *tt;
    tt = &(solid0->IEN[elenum][0]);
      
    for(ii=0;ii<nlbf;ii++)
    {
       row = tt[ii];
       //cout << row << '\t' << rhslocal(ii) << endl;
       rhsVec(row) += rhslocal(ii);
       for(jj=0;jj<nlbf;jj++)
       {
          col = tt[jj];
          
          //cout << row << '\t' << col << '\t' << endl;
          //cout << '\t' << coeffMat.coeff(row,col) << endl;
          //cout << '\t' << coeffMat.coeffRef(row,col) << endl;
          
          coeffMat.coeffRef(row, col) += Nlocal(ii, jj);
       }
    }
*/

  return;
}



void  NurbsElem3DStructMixed2fieldStabilised::AssembleElementMatrix(int index, SparseMatrixXd& mtx, int start1, int start2)
{
    int  ii, jj, nn=0, aa, bb, size2, ind, *tt1, *tt2;

    tt1 = &(solid0->LM[elenum][0]);

    for(ii=0;ii<nsize;ii++)
    {
      aa = tt1[ii];
      if( aa != -1)
      {
        for(jj=0;jj<nsize;jj++)
        {
          bb = tt1[jj];
          if(bb != -1)
          {
            //cout << ii << '\t' << jj << '\t' << aa << '\t' << bb << endl;
            mtx.coeffRef(aa, bb) += Klocal(ii, jj);
          }
        }
        //cout << " ii = " << ii << '\t' << aa << endl;
      }
    }
    //cout << " qqqqqqqqqqqq " << endl;

  return;
}









