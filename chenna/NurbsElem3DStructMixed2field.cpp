#include <Eigen/Dense>

#include <math.h>
#include "Debug.h"
#include "FunctionsElement.h"
#include "MpapTime.h"
#include "NurbsElem3DStructMixed2field.h"
#include "NurbsShapeFunctions.h"
#include <assert.h>
#include "ComputerTime.h"
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
extern PlotVTK plotvtk;


NurbsElem3DStructMixed2field::NurbsElem3DStructMixed2field(void)
{
  if (debug) cout << " constructor NurbsElem3DStructMixed2field\n\n";

//   cout << " constructor NurbsElem3DStructMixed2field\n\n";

  calcExtraMatrices = true;
}



NurbsElem3DStructMixed2field::~NurbsElem3DStructMixed2field()
{
  if (debug) cout << " destructor NurbsElem3DStructMixed2field\n\n";

  // cout << " destructor NurbsElem3DStructMixed2field\n\n";

}

int NurbsElem3DStructMixed2field::calcStiffnessMatrix(double dt)
{

  return 0;
}



void NurbsElem3DStructMixed2field::contourplot(int index, int nCol, double umin, double umax)
{
  return;
}

int NurbsElem3DStructMixed2field::calcMassMatrix(int lumpInd, double dt)
{
  return 0;
}




int NurbsElem3DStructMixed2field::calcStiffnessAndResidual()
{
  if(finite)
    return  NurbsElem3DStructMixed2field::calcStiffnessAndResidual2();
  else
    return  NurbsElem3DStructMixed2field::calcStiffnessAndResidual1();
}




int NurbsElem3DStructMixed2field::calcStiffnessAndResidual1()
{
   int  err, isw, count, count1, ll, ii, jj, kk, index, gp1, gp2, gp3, mm;

   int  TI, TIp1, TIp2, TJ, TJp1, TJp2;
   int  sizep = solid2->nlbf;
   
   double  F[9], detF, fact, dvol0, dt, Jac, pres, bb1, bb2, bb3, utemp, vtemp, wtemp;
   double  eps, volstr, BULK, fact1, fact2;

   double  cc[6][6], stre[6], bc[3][6], Idev[6][6], cctmp[6][6];
   double  N[nlbf], dN_dx[nlbf], dN_dy[nlbf],  dN_dz[nlbf], Nbar[sizep];

   Idev3D(Idev);

   double *gaussweights = &(solid0->gaussweights[0]);

   BULK = matDat[0];

   eps = 1.0/BULK;
   //eps = 0.0;

   //cout << " eps  = " << eps << endl;

   Klocal.setZero();
   Flocal.setZero();
   resi2.zero();

   Kup.resize(nsize, sizep);
   Kup.setZero();

   Kpp.resize(sizep, sizep);
   Kpp.setZero();

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
          
          volstr = (F[0]+F[4]+F[8]-3.0);
          
          //cout  << volstr << '\t' << Nbar[0] << '\t' << dvol0 << endl;

          matlib3d_(matDat, F, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &isw, &err, &count, (Element*) this);
          if(err !=0)           return 1;

          pres = solid2->computeValueAndShanpeFns(1, utemp, vtemp, wtemp, Nbar);

          //cout << " pres " << pres << endl;

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

          for(ii=0;ii<6;ii++)
          {
             stre[ii] *= dvol0;
             for(jj=0;jj<6;jj++)
               cc[ii][jj] *= dvol0;
          }

          //==============================================
          // CALCULATE TANGENT STIFFNESS and RESIDUAL
          //==============================================

            for(ii=0;ii<nlbf;ii++)
            {
                bb1 = dN_dx[ii];
                bb2 = dN_dy[ii];
                bb3 = dN_dz[ii];
                
                for(kk=0;kk<6;kk++)
                {
                   bc[0][kk] = (bb1 * cc[0][kk] + bb2 * cc[3][kk] + bb3 * cc[5][kk]);
                   bc[1][kk] = (bb1 * cc[3][kk] + bb2 * cc[1][kk] + bb3 * cc[4][kk]);
                   bc[2][kk] = (bb1 * cc[5][kk] + bb2 * cc[4][kk] + bb3 * cc[2][kk]);
                }
           
                TI   = 3*ii;
                TIp1 = TI+1;
                TIp2 = TI+2;

                Flocal[TI]   -= (bb1*stre[0] + bb2*stre[3] + bb3*stre[5]) ;
                Flocal[TIp1] -= (bb1*stre[3] + bb2*stre[1] + bb3*stre[4]) ;
                Flocal[TIp2] -= (bb1*stre[5] + bb2*stre[4] + bb3*stre[2]) ;

                for(jj=0;jj<nlbf;jj++)
                {
                   bb1 = dN_dx[jj];
                   bb2 = dN_dy[jj];
                   bb3 = dN_dz[jj];

                   TJ   = 3*jj;
                   TJp1 = TJ+1;
                   TJp2 = TJ+2;

                   Klocal(TI,   TJ)    +=  (bc[0][0] * bb1 + bc[0][3] * bb2 + bc[0][5] * bb3) ;
                   Klocal(TIp1, TJ)    +=  (bc[1][0] * bb1 + bc[1][3] * bb2 + bc[1][5] * bb3) ;
                   Klocal(TIp2, TJ)    +=  (bc[2][0] * bb1 + bc[2][3] * bb2 + bc[2][5] * bb3) ;
                   
                   Klocal(TI,   TJp1)  +=  (bc[0][1] * bb2 + bc[0][3] * bb1 + bc[0][4] * bb3) ;
                   Klocal(TIp1, TJp1)  +=  (bc[1][1] * bb2 + bc[1][3] * bb1 + bc[1][4] * bb3) ;
                   Klocal(TIp2, TJp1)  +=  (bc[2][1] * bb2 + bc[2][3] * bb1 + bc[2][4] * bb3) ;
                   
                   Klocal(TI,   TJp2)  +=  (bc[0][2] * bb3 + bc[0][4] * bb2 + bc[0][5] * bb1) ;
                   Klocal(TIp1, TJp2)  +=  (bc[1][2] * bb3 + bc[1][4] * bb2 + bc[1][5] * bb1) ;
                   Klocal(TIp2, TJp2)  +=  (bc[2][2] * bb3 + bc[2][4] * bb2 + bc[2][5] * bb1) ;
                }

                for(jj=0;jj<sizep;jj++)
                {
                   fact = Nbar[jj] * dvol0;

                   Kup(TI,   jj)  += ( fact * dN_dx[ii] );
                   Kup(TIp1, jj)  += ( fact * dN_dy[ii] );
                   Kup(TIp2, jj)  += ( fact * dN_dz[ii] );
                }
            }

            fact = (volstr - pres*eps)*dvol0;
            fact1 = dvol0*eps;

            for(ii=0;ii<sizep;ii++)
            {
              fact2 = Nbar[ii] * fact1;
              resi2[ii] -= fact * Nbar[ii];

              for(jj=0;jj<sizep;jj++)
                Kpp(ii,jj) -= ( fact2 * Nbar[jj] );
            }

          count++;
          count1++;
          ll += nivGP;
   } // gp1
   } // gp2
   } // gp3

  //printMatrix(Klocal); printf("\n\n");    printVector(Flocal);

  return 0;
}





int NurbsElem3DStructMixed2field::calcStiffnessAndResidual2()
{

   int  err, isw, count, count1, ll, ii, jj, kk, index, gp1, gp2, gp3, mm;
   int  TI, TIp1, TIp2, TJ, TJp1, TJp2;
   int  sizep = solid2->nlbf;

   double  F[9], detF, dvol0, dt, Jac, dummy, pres, bb1, bb2, bb3, utemp, vtemp, wtemp, r1d3 = 1.0/3.0, r2d3 = 2.0*r1d3;
   double  fact, fact1, fact2, fact3, pbar, dvol, BULK, eps;

   double  D11[6][6], stre[6], bc[3][6], Idev[6][6], cctmp[6][6], strdev[6];
   double  N[nlbf], dN_dx[nlbf], dN_dy[nlbf],  dN_dz[nlbf], Nbar[sizep];

   Idev3D(Idev);

   double *gaussweights = &(solid0->gaussweights[0]);

   BULK = matDat[0];

   eps = 1.0/BULK;
   //eps = 0.0;
   
   //cout << " eps  = " << eps << endl;

   Klocal.setZero();
   Flocal.setZero();
   resi2.zero();

   Kup.resize(nsize, sizep);
   Kup.setZero();

   Kpp.resize(sizep, sizep);
   Kpp.setZero();

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

          solid1->ShapeFunDerivatives(&startindex[0], &(knotsAtGPs[index]), N, dN_dx, dN_dy, dN_dz, Jac);

          dvol = Jac * gaussweights[count1] * JacMultFact;

          matlib3d_(matDat, F, stre, D11[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &isw, &err, &count, (Element*) this);
          if(err !=0)           return 1;

          pres = solid2->computeValueAndShanpeFns(1, utemp, vtemp, wtemp, Nbar);

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

          for(ii=0;ii<6;ii++)
          {
             stre[ii] *= dvol;
             for(jj=0;jj<6;jj++)
               D11[ii][jj] *= dvol;
          }

          //==============================================
          // CALCULATE TANGENT STIFFNESS and RESIDUAL
          //==============================================

            for(ii=0;ii<nlbf;ii++)
            {
                bb1 = dN_dx[ii];
                bb2 = dN_dy[ii];
                bb3 = dN_dz[ii];
                
                for(kk=0;kk<6;kk++)
                {
                   bc[0][kk] = (bb1 * D11[0][kk] + bb2 * D11[3][kk] + bb3 * D11[5][kk]);
                   bc[1][kk] = (bb1 * D11[3][kk] + bb2 * D11[1][kk] + bb3 * D11[4][kk]);
                   bc[2][kk] = (bb1 * D11[5][kk] + bb2 * D11[4][kk] + bb3 * D11[2][kk]);
                }

                TI   = 3*ii;
                TIp1 = TI+1;
                TIp2 = TI+2;

                fact1 = bb1 * dvol;
                fact2 = bb2 * dvol;
                fact3 = bb3 * dvol;

                for(jj=0;jj<sizep;jj++)
                {
                   fact = Nbar[jj];

                   Kup(TI,   jj) += ( fact * fact1 );
                   Kup(TIp1, jj) += ( fact * fact2 );
                   Kup(TIp2, jj) += ( fact * fact3 );
                }

                fact1 = bb1*stre[0] + bb2*stre[3] + bb3*stre[5] ;
                fact2 = bb1*stre[3] + bb2*stre[1] + bb3*stre[4] ;
                fact3 = bb1*stre[5] + bb2*stre[4] + bb3*stre[2] ;

                Flocal[TI]   -= fact1 ;
                Flocal[TIp1] -= fact2 ;
                Flocal[TIp2] -= fact3 ;

                for(jj=0;jj<nlbf;jj++)
                {
                   bb1 = dN_dx[jj];
                   bb2 = dN_dy[jj];
                   bb3 = dN_dz[jj];

                   TJ   = 3*jj;
                   TJp1 = TJ+1;
                   TJp2 = TJ+2;

                   fact = fact1 * bb1 + fact2 * bb2 + fact3 * bb3;

                   Klocal(TI,   TJ)    +=  (bc[0][0] * bb1 + bc[0][3] * bb2 + bc[0][5] * bb3 + fact) ;
                   Klocal(TIp1, TJ)    +=  (bc[1][0] * bb1 + bc[1][3] * bb2 + bc[1][5] * bb3) ;
                   Klocal(TIp2, TJ)    +=  (bc[2][0] * bb1 + bc[2][3] * bb2 + bc[2][5] * bb3) ;
                   
                   Klocal(TI,   TJp1)  +=  (bc[0][1] * bb2 + bc[0][3] * bb1 + bc[0][4] * bb3) ;
                   Klocal(TIp1, TJp1)  +=  (bc[1][1] * bb2 + bc[1][3] * bb1 + bc[1][4] * bb3 + fact) ;
                   Klocal(TIp2, TJp1)  +=  (bc[2][1] * bb2 + bc[2][3] * bb1 + bc[2][4] * bb3) ;
                   
                   Klocal(TI,   TJp2)  +=  (bc[0][2] * bb3 + bc[0][4] * bb2 + bc[0][5] * bb1) ;
                   Klocal(TIp1, TJp2)  +=  (bc[1][2] * bb3 + bc[1][4] * bb2 + bc[1][5] * bb1) ;
                   Klocal(TIp2, TJp2)  +=  (bc[2][2] * bb3 + bc[2][4] * bb2 + bc[2][5] * bb1 + fact) ;
                }
            }

            fact1 = dvol0*eps;
            fact  = (detF - 1.0 - pres*eps) * dvol0;

            for(ii=0;ii<sizep;ii++)
            {
               resi2[ii] -= Nbar[ii] * fact;

               fact2 = Nbar[ii] * fact1;

               for(jj=0;jj<sizep;jj++)
                 Kpp(ii,jj) -= ( fact2 * Nbar[jj] );
            }

          count++;
          count1++;
          ll += nivGP;
   } // gp1
   } // gp2
   } // gp3

  return 0;
}





int NurbsElem3DStructMixed2field::calcInternalForces()
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

//  printForceVector();
*/
  return 0;
}




int NurbsElem3DStructMixed2field::calcOutput(double u1, double v1)
{
  return 0;
}




void NurbsElem3DStructMixed2field::discreteContourplot(int vartype, int varindex, int index, int nCol, double umin, double umax)
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

    mapper1->SetInputConnection(uGrid->GetProducerPort());
    //mapper1->SetInputData(uGrid);

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






void NurbsElem3DStructMixed2field::projectToKnots(bool extrapolateFlag, int vartype, int varindex, int index)
{
/*
   vals2project[0] = intVar2[indx];
   vals2project[1] = intVar2[(nGP1-1)*nivGP+indx];
   vals2project[2] = intVar2[nGP1*(nGP2-1)*nivGP+indx];
   vals2project[3] = intVar2[(nGP1*nGP2-1)*nivGP+indx];
*/

   //double outval[nGP];
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




void NurbsElem3DStructMixed2field::projectStress(int varindex, double* outval)
{
/*
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

          pres = solid2->computeValue(1, knotsAtGPs[index], knotsAtGPs[index+1], knotsAtGPs[index+2]);

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
*/
  return;
}




void NurbsElem3DStructMixed2field::projectStrain(int vartype, int varindex, double* outval)
{
  return;
}






void NurbsElem3DStructMixed2field::projectIntVar(int index, double* outval)
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




void NurbsElem3DStructMixed2field::AssembleElementMatrix2(int index, MatrixXd& stiff, MatrixXd& CC)
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



void NurbsElem3DStructMixed2field::AssembleElementVector2(bool firstIter, int flag, VectorXd& rhs, double* reac, VectorXd& uu)
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
      cerr << " ERROR in flag in 'NurbsElem3DStructMixed2field::AssembleElementVector2' " << endl;
      return;
   }


//  cout << " resi " << resi << endl;

  return;
}






void NurbsElem3DStructMixed2field::AssembleElementMatrix(int index, MatrixSparseArray<double>& mtx)
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





void NurbsElem3DStructMixed2field::AssembleElementVector(bool firstIter, bool flag, double* rhs, double* reac, int start1, int start2)
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
      int aa, bb, size2, ind, *tt2;

      tt2   = &(solid2->LM[elenum2][0]);
      size2 = solid2->nsize;

      for(aa=0;aa<nsize;aa++)
      {
         if(tt[aa] != -1)
           rhs[tt[aa]] += Flocal[aa];

         // add up reaction forces
         reac[forassy[aa]] += Flocal[aa];
      }

      for(aa=0;aa<size2;aa++)
      {
         if(tt2[aa] != -1)
           rhs[tt2[aa] + start1] += resi2[aa];
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

         for(aa=0;aa<size2;aa++)
         {
             ind = tt2[aa] + start1;

             for(bb=0;bb<nsize;bb++)
             {
                 if(tt[bb] == -1)
                   rhs[ind] -= Kup(bb,aa) * primvar[bb];
             }
         }
      }
      tt2 = NULL;
   }

//  cout << " resi " << resi << endl;

  tt= NULL;

  return;
}



/*
void  NurbsElem3DStructMixed2field::AssembleElementMatrix(int index, Mat mtx, int start1, int start2)
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



void NurbsElem3DStructMixed2field::toPostprocess(int vartype, int varindex, int type, SparseMatrixXd&  coeffMat, VectorXd& rhsVec)
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



void  NurbsElem3DStructMixed2field::AssembleElementMatrix(int index, SparseMatrixXd& mtx, int start1, int start2)
{
    int  ii, jj, nn=0, aa, bb, size2, ind, *tt1, *tt2;

    tt1 = &(solid0->LM[elenum][0]);
    tt2 = &(solid2->LM[elenum2][0]);
    size2 = solid2->nsize;

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

        for(jj=0;jj<size2;jj++)
        {
          bb = tt2[jj];
          if(bb != -1)
          {
            bb += start1;
            mtx.coeffRef(aa, bb) += Kup(ii, jj);
            mtx.coeffRef(bb, aa) += Kup(ii, jj);
            //mtx.coeffRef(bb, aa) += Kup(ii, jj);
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
            mtx.coeffRef(aa, bb) += Kpp(ii, jj);
          }
        }
      }
    }


  return;
}









