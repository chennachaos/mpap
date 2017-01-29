#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
//#define EIGEN_SUPERLU_SUPPORT
//#define EIGEN_UMFPACK_SUPPORT

#include <iostream>

#include "IsogeometricFEM.h"
#include "FunctionsProgram.h"
#include "DataBlockTemplate.h"
#include "PropertyTypeEnum.h"
#include "MathGeom.h"
#include "NurbsShapeFns.h"
#include "NurbsMiscFunctions.h"
#include "ComputerTime.h"
#include "MpapTime.h"
//#include "Plot.h"
#include <assert.h>
#include "PlotVTK.h"
#include "SolverMA41Eigen.h"

#include <fstream>
#include <string.h>
#include <iomanip>
#include <iostream>
#include <fstream>


#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
//#include <Eigen/UmfPackSupport>
//#include <Eigen/SuperLUSupport>
#include <Eigen/IterativeSolvers>

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
#include <vtkScalarBarActor.h>
#include "vtkActor.h"
#include "vtkConeSource.h"
#include "vtkGlyph3D.h"
#include <vtkLookupTable.h>
#include "vtkPolyDataMapper.h"
#include "vtkRenderer.h"
#include "vtkSphereSource.h"
#include "vtkXOpenGLRenderWindow.h"
#include "vtkXRenderWindowInteractor.h"
#include "vtkTextProperty.h"
#include "vtkContourFilter.h"
#include <vtkQuad.h>

#include "Files.h"

//#include <Xm/PushB.h>
//#include <Xm/Form.h>

//extern Plot plot;
extern PlotVTK  plotvtk;
extern ComputerTime computerTime;
extern MpapTime     mpapTime;
extern Files     files;


using namespace std;
using namespace Eigen;


//typedef SparseMatrix<double> SparseMatrixXd;



void IsogeometricFEM::CalcConstraintVariables(int index)
{
  int ind1=0, n1, n2, ii, jj, iii;

  for(iii=0;iii<Npatch;iii++)
  {
    n1 = uu[iii].n;
    n2 = vv[iii].n;
    for(jj=0;jj<n2;jj++)
    {
      ind1 = n1*jj;
      for(ii=0;ii<n1;ii++)
      {
        Sdef[iii][ind1+ii] = SurfaceResult[iii].SurfacePoint(uu[iii][ii], vv[iii][jj]).CalcEuclid();
        outp[iii][ind1+ii] = surfSecondVar[iii].computeValue(1, uu[iii][ii], vv[iii][jj]) ;
      }
    }
  }
  
  //cout << " llllllllllll  " << endl;

  return;
}


void IsogeometricFEM::CalcDisplacements(int dir)
{
  if(dir > 3)
  {
    cout << endl;
    cout << "         This item does not exist for Displacement plot " << endl;
    return;
  }

  assert(RSYS <= 1);

  EPOINT *EP1, *EP2;
  int ind1=0, n1, n2, ii, jj, iii;

  if(RSYS == 0)
  {
    for(iii=0;iii<Npatch;iii++)
    {
      n1 = uu[iii].n;
      n2 = vv[iii].n;
      for(jj=0;jj<n2;jj++)
      {
         ind1 = n1*jj;
         for(ii=0;ii<n1;ii++)
         {
            Sdef[iii][ind1+ii] = SurfaceResult[iii].SurfacePoint(uu[iii][ii], vv[iii][jj]).CalcEuclid();
            //Sdef[iii][ind1+ii] = SurfaceResult[iii].SurfacePoint2(uu[iii][ii], vv[iii][jj]);

            EP1 = &(Sorig[iii][ind1+ii]);
            EP2 = &(Sdef[iii][ind1+ii]);

            if(dir == 0)
              outp[iii][ind1+ii] = EP2->x - EP1->x ;
            else if(dir == 1)
              outp[iii][ind1+ii] = EP2->y - EP1->y ;
            else
              outp[iii][ind1+ii] = EP2->z;
         }
      }
    }
  }
  else
  {
    double f1, f2;
    for(iii=0;iii<Npatch;iii++)
    {
      n1 = uu[iii].n;
      n2 = vv[iii].n;
      for(jj=0;jj<n2;jj++)
      {
         ind1 = n1*jj;
         for(ii=0;ii<n1;ii++)
         {
            Sdef[iii][ind1+ii] = SurfaceResult[iii].SurfacePoint(uu[iii][ii], vv[iii][jj]).CalcEuclid();
            //Sdef[iii][ind1+ii] = SurfaceResult[iii].SurfacePoint2(uu[iii][ii], vv[iii][jj]);

            EP1 = &(Sorig[iii][ind1+ii]);
            EP2 = &(Sdef[iii][ind1+ii]);

            f1 = EP2->x - EP1->x;
            f2 = EP2->y - EP1->y;

            if(dir == 0)
              outp[iii][ind1+ii] = sqrt(f1*f1 + f2*f2) ;
            else
              outp[iii][ind1+ii] = atan(EP2->y/EP2->x) - atan(EP1->y/EP1->x);
         }
      }
    }
  }

  return;
}



void IsogeometricFEM::contourplotStress(int vartype, int varindex, bool extrapolateFlag, int index, int nCol, bool umnxflag, double umin, double umax, bool legdflag)
{
  cout << "     ISOGEOMETRICFEM: contour plot ... " << endl;  cout << endl;
  return;
}


void IsogeometricFEM::contourplot(int vartype, int varindex, bool extrapolateFlag, int index, int nCol, bool umnxflag, double umin, double umax, bool legdflag)
{
  cout << "     ISOGEOMETRICFEM: contour plot ... " << endl;  cout << endl;

//  SurfaceResult[0].PlotValues(1);
  //return;

  if(patchGrp[0].ndom == 1)
  {
     cout << "        ERROR...: Contour Plot is not available for Curves ... " << endl;
     cout << endl;
     return;
  }
  
  printf("vartype = %5d \t varindex = %5d \t index = %5d \n", vartype, varindex, index);


  switch(vartype)
  {
    case 0:  // plot displacement

          CalcDisplacements(varindex);

          break;

    case 1:  // plot total strain
    case 2:  // plot elastic strain
    case 3:  // plot plastic strain
    case 4:  // plot stress

            if(mixedSolverFlag > 6 && varindex == 5)
              CalcConstraintVariables(index+1);
            else
              projectFromElemsToKnots(extrapolateFlag, vartype, varindex, index);

          break;

    case 5:  // plot element internal variables

             projectFromElemsToKnots(extrapolateFlag, vartype, varindex, index);
          break;

    default:
             cout << endl;
             cout << "         This item does not exist for contour plot " << endl;
             cout << endl;
             break;
  }
/*
 if(vartype > 0)
 {
    SparseMatrixXd  coeffMat(ntotgbf,ntotgbf);
    VectorXd  rhsVec(ntotgbf), soln(ntotgbf);

    coeffMat.setZero();
    rhsVec.setZero();

    for(int ee=0;ee<totnumel;ee++)
      elem[ee]->toPostprocess(vartype, varindex, index, coeffMat, rhsVec);

    //soln = coeffMat.ldlt().solve(rhsVec);
    
    ConjugateGradient<SparseMatrix<double>, Eigen::Upper> solver;
    soln = solver.compute(coeffMat).solve(rhsVec);


    VectorArray<double> tmpVec;  tmpVec.setDim(ntotgbf);
    
    for(int ii=0;ii<ntotgbf;ii++)
      tmpVec[ii] = soln(ii);
      
    SurfaceResult[0].Values[0] = tmpVec;

    CalcConstraintVariables(1);
  }
*/

    double umn = 0.0, umx = 0.0;

    umn = outp[0][0];
    umx = outp[0][0];
    for(int iii=0;iii<Npatch;iii++)
    {
       for(int jj=0;jj<outp[iii].n;jj++)
       {
          //cout << outp[iii][jj] << endl;
          if(outp[iii][jj] < umn) umn = outp[iii][jj];
          if(outp[iii][jj] > umx) umx = outp[iii][jj];
       }
    }

//  if(mixedSolverFlag > 6)
//    cout << surfSecondVar[0].Values << endl; cout << endl;
//    cout << '\t' << outp[0] << endl;

    if(!umnxflag)
    {
      umin = umn;
      umax = umx;
    }

//    cout << '\t' << " umin and umax  " << umin << '\t' << umax << endl;

     //////////////////////////////////////////
     //                                      //
     //          3---4                       //
     //          |  /|                       //
     //          | / |                       //
     //          |/  |                       //
     //          1---2                       //
     //                                      //
     //////////////////////////////////////////

    EPOINT *EP1, *EP2, *EP3, *EP4;
    double x1[2], x2[2], x3[2], x4[2];
    double u1, u2, u3, u4;
    int ind1=0, ind2=0, n1, n2, iii, ii, jj, temp;

    bool  finite = (patchElemProp[0].data[2] == 1);

    //cout << finite << endl;

    if(defUndefFlag)
    {
        for(iii=0;iii<Npatch;iii++)
        {
            n1 = uu[iii].n;
            n2 = vv[iii].n;

            for(jj=0;jj<(n2-1);jj++)
            {
               ind1 = n1*jj;
               ind2 = n1*(jj+1);
               for(ii=0;ii<(n1-1);ii++)
               {
                  temp = ind1+ii;
                  EP1 = &(Sdef[iii][temp]) ;
                  u1 = outp[iii][temp];

                  temp += 1;
                  EP2 = &(Sdef[iii][temp]) ;
                  u2 = outp[iii][temp];

                  temp = ind2+ii;
                  EP3 = &(Sdef[iii][temp]) ;
                  u3 = outp[iii][temp];

                  temp += 1;
                  EP4 = &(Sdef[iii][temp]) ;
                  u4 = outp[iii][temp];

                  x1[0] = EP1->x; x1[1] = EP1->y;
                  x2[0] = EP2->x; x2[1] = EP2->y;
                  x3[0] = EP3->x; x3[1] = EP3->y;
                  x4[0] = EP4->x; x4[1] = EP4->y;

                  // contour plot for 1st triangle
                  ////plot.triangleContourPlot(x1, x2, x4, u1, u2, u4, umin, umax, nCol);

                  // contour plot for 2nd triangle
                  ////plot.triangleContourPlot(x1, x3, x4, u1, u3, u4, umin, umax, nCol);
               }
	    }
        }
    }
    else
    {
        for(iii=0;iii<Npatch;iii++)
        {
            n1 = uu[iii].n;
            n2 = vv[iii].n;

            for(jj=0;jj<(n2-1);jj++)
            {
               ind1 = n1*jj;
               ind2 = n1*(jj+1);
               for(ii=0;ii<(n1-1);ii++)
               {
                  temp = ind1+ii;
                  EP1 = &(Sorig[iii][temp]) ;
                  u1 = outp[iii][temp];

                  temp += 1;
                  EP2 = &(Sorig[iii][temp]) ;
                  u2 = outp[iii][temp];

                  temp = ind2+ii;
                  EP3 = &(Sorig[iii][temp]) ;
                  u3 = outp[iii][temp];

                  temp += 1;
                  EP4 = &(Sorig[iii][temp]) ;
                  u4 = outp[iii][temp];

                  x1[0] = EP1->x; x1[1] = EP1->y;
                  x2[0] = EP2->x; x2[1] = EP2->y;
                  x3[0] = EP3->x; x3[1] = EP3->y;
                  x4[0] = EP4->x; x4[1] = EP4->y;

                  // contour plot for 1st triangle
                  ////plot.triangleContourPlot(x1, x2, x4, u1, u2, u4, umin, umax, nCol);

                  // contour plot for 2nd triangle
                  ////plot.triangleContourPlot(x1, x4, x3, u1, u4, u3, umin, umax, nCol);
               }
	    }
        }
    }

    //plot.contourPlotLegend(umin,umax,umn,umx,umnxflag,nCol);

  return;
}



void IsogeometricFEM::projectFromElemsToKnots(bool extrapolateFlag, int vartype, int varindex, int index)
{
   if(patchGrp[0].ndom == 2)
   {
      for(int ee=0;ee<totnumel;ee++)
         elem[ee]->projectToKnots(extrapolateFlag, vartype, varindex, index);

       int nel1, nel2, n1, n2, ind, ind1, ind2, ind3, ind4;
       int iii, e1, e2, ee, ii, jj;

       // zero 'outp'
       for(iii=0;iii<patchGrp.n;iii++)
         outp[iii].zero();

       int cnt = 0;
       ee = 0;
       for(iii=0;iii<patch.n;iii++)
       {
          nel1 = SurfaceListFinal[iii].nelem1;
          nel2 = SurfaceListFinal[iii].nelem2;
          n1 = nel1+1;
          n2 = nel2+1;

          for(e2=0;e2<nel2;e2++)
          {
             //ind = nel1*e2;
             for(e1=0;e1<nel1;e1++)
             {
                 ind1 = n1*e2+e1; ind2 = ind1+1;
                 ind3 = n1*(e2+1)+e1; ind4 = ind3+1;
                 //ee = ind+e1;

                 outp[iii][ind1] += elem[ee]->vals2project[0];
                 outp[iii][ind2] += elem[ee]->vals2project[1];
                 outp[iii][ind3] += elem[ee]->vals2project[2];
                 outp[iii][ind4] += elem[ee]->vals2project[3];

                 ee++;
             }
          }
       }

        for(iii=0;iii<patch.n;iii++)
        {
            n1 = uu[iii].n;
            n2 = vv[iii].n;

            for(jj=0;jj<n2;jj++)
            {
                ind1 = n1*jj;
                if(jj==0)
                {
                    for(ii=1;ii<(n1-1);ii++)
                    {
                       outp[iii][ind1+ii] /= 2.0;
                    }
                }
                else if( jj==(n2-1) )
                {
                    for(ii=1;ii<(n1-1);ii++)
                    {
                       outp[iii][ind1+ii] /= 2.0;
                    }
                }
                else
                {
                    outp[iii][ind1] /= 2.0; //ii=0
                    outp[iii][ind1+n1-1] /= 2.0; // ii= n1-1
                    for(ii=1;ii<(n1-1);ii++)
                    {
                       outp[iii][ind1+ii] /= 4.0;
                    }
                }
            }
        }

        //
        if(patch.n > 1)
        {
            double val1, val2;
            int ind2, ee1, ee2, ee3, ee4, nelm11, nelm12, nelm21, nelm22, nn1=0, nn2=0, jj;

            nn1 = 0;
            nn2 = SurfaceListFinal[0].nelem;

            for(iii=0;iii<patch.n-1;iii++)
            {
                n1 = uu[iii].n;
                n2 = uu[iii+1].n;

                nelm11 = SurfaceListFinal[iii].nelem1;
                nelm12 = SurfaceListFinal[iii].nelem2;
                
                nelm21 = SurfaceListFinal[iii+1].nelem1;
                nelm22 = SurfaceListFinal[iii+1].nelem2;
                
                nn1 = (nelm11+1)*nelm12;
                nn2 = 0;

                for(jj=0; jj<=nelm11; jj++)
                {
                  ind1 = nn1+jj;
                  ind2 = nn2+jj;

                  val1 = 0.5*(outp[iii][ind1] + outp[iii+1][ind2] );

                  outp[iii][ind1] = val1;
                  outp[iii+1][ind2] = val1;
                }

                nn1 += SurfaceListFinal[iii].nelem;
                nn2 += SurfaceListFinal[iii+1].nelem;
            }
        }
        //
        /*
        if(patch.n > 1)
        {
            double val1, val2;
            int ind2, ee1, ee2, ee3, ee4, nelm1, nelm2, nn1=0, nn2=0, jj;

            nn1 = 0;
            nn2 = SurfaceListFinal[0].nelem;

            for(iii=0;iii<patch.n-1;iii++)
            {
                n1 = uu[iii].n;
                n2 = uu[iii+1].n;

                nelm1 = SurfaceListFinal[iii].nelem1;
                nelm2 = SurfaceListFinal[iii+1].nelem1;

                // for the case jj = 0
                jj = 0;
                ind1 = n1*(jj+1) - 1;
                ind2 = n2*jj;

                ee1 = nn1 + nelm1 - 1;
                ee2 = nn2;

                val1 = 0.5*(elem[ee1]->vals2project[1] + elem[ee2]->vals2project[0]);

                outp[iii][ind1] = val1;
                outp[iii+1][ind2] = val1;

	         // for the case jj = vv[iii].n-1

	         jj = vv[iii].n-1;
                ind1 = n1*(jj+1) - 1;
                ind2 = n2*jj;

                ee1 = nn1 + nelm1*jj - 1;
                ee2 = nn2 + nelm2*(jj-1) ;

                val1 = 0.5*(elem[ee1]->vals2project[3] + elem[ee2]->vals2project[2]);

                outp[iii][ind1] = val1;
                outp[iii+1][ind2] = val1;

                for(jj=1;jj<vv[iii].n-1;jj++)
                {
                    ind1 = n1*(jj+1) - 1;
                    ind2 = n2*jj;

                    ee1 = nn1 + (nelm1*jj) - 1;
                    ee2 = nn1 + (nelm1*(jj+1)) - 1;
                    ee3 = nn2 + (nelm2*(jj-1)) ;
                    ee4 = nn2 + (nelm2*jj) ;

                    val1 = 0.25*(elem[ee1]->vals2project[3] + elem[ee2]->vals2project[1] + elem[ee3]->vals2project[2] + elem[ee4]->vals2project[0]);

               //     cout << '\t' << jj << '\t' << elem[ee1]->vals2project[3] << '\t' << elem[ee2]->vals2project[1] << '\t' << elem[ee3]->vals2project[2] << '\t' << elem[ee4]->vals2project[0] << endl;

                    outp[iii][ind1] = val1;
                    outp[iii+1][ind2] = val1;
                }

                nn1 += SurfaceListFinal[iii].nelem;
                nn2 += SurfaceListFinal[iii+1].nelem;
            }
        }
        */
     }

   if(patchGrp[0].ndom == 3)
   {
       int nel1, nel2, nel3, n1, n2, n3, ind, ind1, ind2, ind3, ind4, nn;
       int iii, e1, e2, e3, ee, ii, jj, kk, plane1, plane2, ind5, ind6, ind7, ind8;
   
       double *tt;

       for(ee=0;ee<totnumel;ee++)
          elem[ee]->projectToKnots(extrapolateFlag, vartype, varindex, index);

       // zero 'outp'
       for(iii=0;iii<patchGrp.n;iii++)
         outp[iii].zero();

       int cnt = 0;
       ee = 0;
       for(iii=0;iii<patch.n;iii++)
       {
          nel1 = SolidListFinal[iii].nelem1;
          nel2 = SolidListFinal[iii].nelem2;
          nel3 = SolidListFinal[iii].nelem3;

          n1 = nel1+1;      n2 = nel2+1;      n3 = nel3+1;
      
          nn = n1*n2;

          for(e3=0;e3<nel3;e3++)
          {
              plane1 = nn*e3;
              plane2 = nn*(e3+1);
              for(e2=0;e2<nel2;e2++)
              {
                 //ind = nel1*e2;
                 for(e1=0;e1<nel1;e1++)
                 {
                     ind1 = plane1 + n1*e2+e1;      ind2 = ind1+1;
                     ind3 = plane1 + n1*(e2+1)+e1;  ind4 = ind3+1;

                     ind5 = plane2 + n1*e2+e1;      ind6 = ind5+1;
                     ind7 = plane2 + n1*(e2+1)+e1;  ind8 = ind7+1;
                 
                     tt = &(elem[ee]->vals2project[0]);

                     outp[iii][ind1] += tt[0];
                     outp[iii][ind2] += tt[1];
                     outp[iii][ind3] += tt[2];
                     outp[iii][ind4] += tt[3];

                     outp[iii][ind5] += tt[4];
                     outp[iii][ind6] += tt[5];
                     outp[iii][ind7] += tt[6];
                     outp[iii][ind8] += tt[7];

                     ee++;
                 }
              }
          }
       }

        for(iii=0;iii<patch.n;iii++)
        {
            n1 = uu[iii].n;
            n2 = vv[iii].n;
            n3 = ww[iii].n;
        
            nn = n1*n2;
        
            for(kk=0;kk<n3;kk++)
            {
                plane1 = nn*kk;
                if( kk==0 || kk==(n3-1) )
                {
                    for(jj=0;jj<n2;jj++)
                    {
                        ind1 = plane1+n1*jj;
                        if(jj==0 || jj==(n2-1) )
                        {
                            for(ii=1;ii<(n1-1);ii++)
                            {
                               outp[iii][ind1+ii] /= 2.0;
                            }
                        }
                        else
                        {
                            outp[iii][ind1] /= 2.0; //ii=0
                            outp[iii][ind1+n1-1] /= 2.0; // ii= n1-1
                            for(ii=1;ii<(n1-1);ii++)
                            {
                               outp[iii][ind1+ii] /= 4.0;
                            }
                        }
                    }
                }
                else
                {
                    for(jj=0;jj<n2;jj++)
                    {
                        ind1 = plane1+n1*jj;
                        if(jj==0 || jj==(n2-1) )
                        {
                            outp[iii][ind1] /= 2.0; //ii=0
                            outp[iii][ind1+n1-1] /= 2.0; // ii= n1-1
                            for(ii=1;ii<(n1-1);ii++)
                            {
                               outp[iii][ind1+ii] /= 4.0;
                            }
                        }
                        else
                        {
                            outp[iii][ind1] /= 4.0; //ii=0
                            outp[iii][ind1+n1-1] /= 4.0; // ii= n1-1
                            for(ii=1;ii<(n1-1);ii++)
                            {
                               outp[iii][ind1+ii] /= 8.0;
                            }
                        }
                    }
                }
            }
        }
   }

  return;
}




void IsogeometricFEM::dispContourplotVTK(int index, int nCol, bool umnxflag, double umin, double umax, int* resln)
{
    plotvtk.clearWindow();

    int   iii, ii, jj, kk, ll, ind1, ind2, ind3, ind4, ind5, ind6;
    vector<int> n1(Npatch), n2(Npatch), n3(Npatch), nn(Npatch);

    ListArray<VectorArray<double> >  uu, vv, ww;

    vtkIdType  N1, N2, count;
    
    uu.setDim(Npatch);
    vv.setDim(Npatch);
    ww.setDim(Npatch);

    count = 0;
    for(iii=0;iii<Npatch;iii++)
    {
       create_vector2(SolidListFinal[iii].U, resln[0], uu[iii]);
       create_vector2(SolidListFinal[iii].V, resln[1], vv[iii]);
       create_vector2(SolidListFinal[iii].W, resln[2], ww[iii]);
       n1[iii] = uu[iii].n;
       n2[iii] = vv[iii].n;
       n3[iii] = ww[iii].n;
       
       nn[iii] = n1[iii] * n2[iii] * n3[iii];
       count += nn[iii];
    }

//    cout << resln[0] << '\t' << resln[1] << '\t' << resln[2] << endl;

    vtkSmartPointer<vtkPoints>            points    =   vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkFloatArray>        scalars   =   vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkScalarBarActor>    scalarBar =   vtkSmartPointer<vtkScalarBarActor>::New();
    vtkSmartPointer<vtkLookupTable>       hueLut    =   vtkSmartPointer<vtkLookupTable>::New();
    vtkSmartPointer<vtkDataSetMapper>     mapper1   =   vtkSmartPointer<vtkDataSetMapper>::New();
    vtkSmartPointer<vtkActor>             actor1    =   vtkSmartPointer<vtkActor>::New();
    vtkSmartPointer<vtkHexahedron>        hex       =   vtkSmartPointer<vtkHexahedron>::New();

    static vtkIdType pts[8];

    EPOINT  EP1, EP2;

    scalars->SetNumberOfTuples(count);
    
    float  val;
    double umn = 0.0, umx = 0.0;
    
    if(index == 0)
    {
        scalars->SetName("UX");
        scalarBar->SetTitle("UX");
    }
    else if(index == 1)
    {
        scalars->SetName("UY");
        scalarBar->SetTitle("UY");
    }
    else
    {
        scalars->SetName("UZ");
        scalarBar->SetTitle("UZ");
    }

    count=0; 
    for(iii=0;iii<Npatch;iii++)
    {    
         for(kk=0;kk<n3[iii];kk++)
         {
             for(jj=0;jj<n2[iii];jj++)
             {
                 for(ii=0;ii<n1[iii];ii++)
                 {
                     EP1 = SolidListFinal[iii].SolidPoint(uu[iii][ii], vv[iii][jj], ww[iii][kk]).CalcEuclid();
                     EP2 = SolidResult[iii].SolidPoint(uu[iii][ii], vv[iii][jj], ww[iii][kk]).CalcEuclid();
              
                     points->InsertNextPoint(EP2.x, EP2.y, EP2.z);
              
                     if(index == 0)
                       val = EP2.x-EP1.x;
                     else  if(index == 1)
                       val = EP2.y-EP1.y;
                     else // if(index == 2)
                       val = EP2.z-EP1.z;

                     if(val < umn) umn = val;
                     if(val > umx) umx = val;
          
                     scalars->SetTuple1(count++, val);
                  }
              }
         }
    }

    if(umnxflag)
    {
      umn = umin;
      umx = umax;
    }
    cout << " umin & umax ... : " << umn << '\t' << umx << endl;
    
    count = 0;
    for(iii=0;iii<Npatch;iii++)
    {    
        for(kk=0;kk<n3[iii]-1;kk++)
        {
            ind5 = count + n1[iii]*n2[iii]*kk;
            ind6 = count + n1[iii]*n2[iii]*(kk+1);

            for(jj=0;jj<n2[iii]-1;jj++)
            {
               ind1 = ind5 + n1[iii]*jj;
               ind2 = ind5 + n1[iii]*(jj+1);
       
               ind3 = ind6 + n1[iii]*jj;
               ind4 = ind6 + n1[iii]*(jj+1);

               for(ii=0;ii<n1[iii]-1;ii++)
               {
                  pts[0] = ind1+ii;          pts[4] = ind3+ii;

                  pts[1] = pts[0]+1;         pts[5] = pts[4]+1;
                  pts[3] = ind2+ii;          pts[7] = ind4+ii;
                  pts[2] = pts[3]+1;         pts[6] = pts[7]+1;
              
                  for(ll=0;ll<8;ll++)
                     hex->GetPointIds()->SetId(ll, pts[ll]);
          
                  plotvtk.uGrid->InsertNextCell(hex->GetCellType(), hex->GetPointIds());
               }
            }
        }
        count += nn[iii];
    }

    mapper1->ScalarVisibilityOn();
    mapper1->SetScalarModeToUsePointData();
    mapper1->SetColorModeToMapScalars();
    mapper1->SetScalarRange(umn, umx);
    mapper1->InterpolateScalarsBeforeMappingOn();

    hueLut->SetHueRange(0.66667, 0.0);
    hueLut->SetTableRange(umn, umx);
    hueLut->SetNumberOfColors(nCol);
    hueLut->SetRampToLinear();
    hueLut->SetScaleToLinear();
    hueLut->Build();
    
    vtkSmartPointer<vtkTextProperty>    txtprop   =   vtkSmartPointer<vtkTextProperty>::New();
    
    txtprop->SetColor(0.0,0.0,0.0);
    txtprop->BoldOn();
    txtprop->ItalicOff();


    mapper1->SetLookupTable(hueLut);

    scalarBar->SetLookupTable(hueLut);
    scalarBar->SetMaximumNumberOfColors(nCol);
    scalarBar->SetNumberOfLabels(nCol+1);
    scalarBar->SetTitleTextProperty(txtprop);
    scalarBar->SetLabelTextProperty(txtprop);
    scalarBar->SetWidth(0.12);
    scalarBar->SetLabelFormat("%10.4f");
    
    plotvtk.uGrid->SetPoints(points);
    plotvtk.uGrid->GetPointData()->SetScalars(scalars);

#if VTK_MAJOR_VERSION == 5
    mapper1->SetInputConnection(plotvtk.uGrid->GetProducerPort());
#else
    mapper1->SetInputData(plotvtk.uGrid);
#endif

    actor1->SetMapper(mapper1);

    plotvtk.rendr->AddActor(actor1);
    plotvtk.rendr->AddActor2D(scalarBar);
    plotvtk.rendr->ResetCamera();
    plotvtk.renWindow->Render();

    return;
}



void IsogeometricFEM::contourplotVTK(int vartype, int varindex, bool extrapolateFlag, int index, int nCol, bool umnxflag, double umin, double umax, int* resln)
{
    plotvtk.clearWindow();

       switch(vartype)
       {
           case 0:  // plot displacement

                   dispContourplotVTK(varindex, nCol, umnxflag, umin, umax, resln);
                   
                   return;

                   break;

           case 1:  // plot total strain
           case 2:  // plot elastic strain
           case 3:  // plot plastic strain
           case 4:  // plot stress
           
                 //  if(mixedSolverFlag > 6 && varindex == 5)
                 //       CalcConstraintVariables(index+1);
                 //  else
                      projectFromElemsToKnots(extrapolateFlag, vartype, varindex, index);
                      
                      //return;

                   break;

           case 5:  // plot element internal variables

                     projectFromElemsToKnots(extrapolateFlag, vartype, varindex, index);
                      //return;
                    break;

           default:
                  cout << endl;
                  cout << "         This item does not exist for contour plot " << endl;
                  cout << endl;
                  break;
       }



    int   ii, jj, kk, ll, n1, n2, n3, ind1, ind2, ind3, ind4, ind5, ind6, nn;

    VectorArray<double>  uu, vv, ww;

    findunique(SolidListFinal[0].U, uu);
    findunique(SolidListFinal[0].V, vv);
    findunique(SolidListFinal[0].W, ww);

    cout << resln[0] << '\t' << resln[1] << '\t' << resln[2] << endl;

    vtkSmartPointer<vtkPoints>            points    =   vtkSmartPointer<vtkPoints>::New();
//    vtkSmartPointer<vtkUnstructuredGrid>  uGrid     =   vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkFloatArray>        scalars   =   vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkScalarBarActor>    scalarBar =   vtkSmartPointer<vtkScalarBarActor>::New();
    vtkSmartPointer<vtkLookupTable>       hueLut    =   vtkSmartPointer<vtkLookupTable>::New();
    vtkSmartPointer<vtkDataSetMapper>     mapper1   =   vtkSmartPointer<vtkDataSetMapper>::New();
    vtkSmartPointer<vtkActor>             actor1    =   vtkSmartPointer<vtkActor>::New();
    vtkSmartPointer<vtkHexahedron>        hex       =   vtkSmartPointer<vtkHexahedron>::New();

    static vtkIdType pts[8];

    EPOINT  EP;

    n1 = uu.n;
    n2 = vv.n;
    n3 = ww.n;
    
    nn = n1*n2;

    scalars->SetNumberOfTuples(nn*n3);
    
    float  val;
    double umn, umx;
    bool  finite = SolidListFinal[0].ElemProp.data[3];
    finite = 1;
    cout << finite << endl;
    
    vtkIdType  N1, N2, count;
    
    umn = umx = outp[0][0];
    
    count=0;
    for(kk=0;kk<n3;kk++)
    {
        for(jj=0;jj<n2;jj++)
        {
            for(ii=0;ii<n1;ii++)
            {
                if(finite)
                   EP = SolidResult[0].SolidPoint(uu[ii], vv[jj], ww[kk]).CalcEuclid();
                else
                   EP = SolidListFinal[0].SolidPoint(uu[ii], vv[jj], ww[kk]).CalcEuclid();
              
                points->InsertNextPoint(EP.x, EP.y, EP.z);
              
                val = outp[0][count];

                if(val < umn) umn = val;
                if(val > umx) umx = val;
          
                scalars->SetTuple1(count, val);
                count++;
             }
         }
    }

    if(umnxflag)
    {
      umn = umin;
      umx = umax;
    }
    cout << " umin & umax ... : " << umn << '\t' << umx << endl;
    
    for(kk=0;kk<n3-1;kk++)
    {
        ind5 = nn*kk;
        ind6 = nn*(kk+1);

        for(jj=0;jj<n2-1;jj++)
        {
           ind1 = ind5 + n1*jj;
           ind2 = ind5 + n1*(jj+1);
       
           ind3 = ind6 + n1*jj;
           ind4 = ind6 + n1*(jj+1);

           for(ii=0;ii<n1-1;ii++)
           {
              pts[0] = ind1+ii;          pts[4] = ind3+ii;

              pts[1] = pts[0]+1;         pts[5] = pts[4]+1;
              pts[3] = ind2+ii;          pts[7] = ind4+ii;
              pts[2] = pts[3]+1;         pts[6] = pts[7]+1;
              
              for(ll=0;ll<8;ll++)
                 hex->GetPointIds()->SetId(ll, pts[ll]);
          
              plotvtk.uGrid->InsertNextCell(hex->GetCellType(), hex->GetPointIds());
           }
        }
    }
    
    plotvtk.uGrid->SetPoints(points);

    plotvtk.uGrid->GetPointData()->SetScalars(scalars);

#if VTK_MAJOR_VERSION == 5
    mapper1->SetInputConnection(plotvtk.uGrid->GetProducerPort());
#else
    mapper1->SetInputData(plotvtk.uGrid);
#endif

//    mapper1->ScalarVisibilityOn();
//    mapper1->SetScalarModeToUsePointData();
//    mapper1->SetColorModeToMapScalars();
    
    mapper1->SetScalarRange(umn, umx);
    actor1->SetMapper(mapper1);
    
    hueLut->SetHueRange(0.66667, 0.0);
    //hueLut->SetTableRange (0, 3);
    //hueLut->SetNumberOfTableValues(nCol);
    hueLut->SetRampToLinear();
    //hueLut->SetScaleToLinear();
    hueLut->Build();
 
    mapper1->SetLookupTable(hueLut);
    scalarBar->SetLookupTable(hueLut);
    scalars->SetLookupTable(hueLut);

    scalarBar->SetMaximumNumberOfColors(nCol);
    scalarBar->SetNumberOfLabels(nCol+1);
    
    scalarBar->SetWidth(0.12);
    scalarBar->SetLabelFormat("%10.4f");
    
    plotvtk.rendr->AddActor(actor1);
    plotvtk.rendr->AddActor2D(scalarBar);
    plotvtk.rendr->ResetCamera();
    plotvtk.renWindow->Render();

  return;
}




void IsogeometricFEM::discreteContourplot(int vartype, int varindex, int index, int nCol, bool umnxflag, double umin, double umax, bool legdflag)
{
   assert(totnumel > 0);

  cout << "     ISOGEOMETRICFEM: discrete contour plot ...\n\n";
  cout << endl;

   for(int ee=0;ee<totnumel;ee++)
   {
    //  cout << '\t' << "     Elem :" << ee << endl;
      elem[ee]->discreteContourplot(vartype, varindex, index, nCol, umin, umax);
   }

   //if(legdflag)
      //plot.contourPlotLegend(umin,umax,umin,umax,umnxflag,nCol);

  return;
}




void IsogeometricFEM::contourplotVTK2(int vartype, int varindex, bool extrapolateFlag, int index, int nCol, bool umnxflag, double umin, double umax, int* resln)
{
       if(vartype > 6 )
       {
          cout << "         This item does not exist for contour plot " << endl;
       }

       plotvtk.clearWindow();

       if(vartype == 0)
       {
           dispContourplotVTK(varindex, nCol, umnxflag, umin, umax, resln);
           return;
       }


      int   iii, ii, jj, kk, ll, ind1, ind2, ind3, ind4, ind5, ind6;
      vector<int>  n1(Npatch), n2(Npatch), n3(Npatch), nn(Npatch);

      ListArray<VectorArray<double> >  uu, vv, ww;

      vtkIdType  count, pts[8];
    
      uu.setDim(Npatch);
      vv.setDim(Npatch);
      ww.setDim(Npatch);

      count = 0;
      for(iii=0;iii<Npatch;iii++)
      {
         create_vector2(SolidListFinal[iii].U, resln[0], uu[iii]);
         create_vector2(SolidListFinal[iii].V, resln[1], vv[iii]);
         create_vector2(SolidListFinal[iii].W, resln[2], ww[iii]);
         n1[iii] = uu[iii].n;
         n2[iii] = vv[iii].n;
         n3[iii] = ww[iii].n;
       
         nn[iii] = n1[iii] * n2[iii] * n3[iii];
         count += nn[iii];
      }

      vtkSmartPointer<vtkPoints>            points    =   vtkSmartPointer<vtkPoints>::New();
      vtkSmartPointer<vtkFloatArray>        scalars   =   vtkSmartPointer<vtkFloatArray>::New();
      vtkSmartPointer<vtkScalarBarActor>    scalarBar =   vtkSmartPointer<vtkScalarBarActor>::New();
      vtkSmartPointer<vtkLookupTable>       hueLut    =   vtkSmartPointer<vtkLookupTable>::New();
      vtkSmartPointer<vtkDataSetMapper>     mapper1   =   vtkSmartPointer<vtkDataSetMapper>::New();
      vtkSmartPointer<vtkActor>             actor1    =   vtkSmartPointer<vtkActor>::New();
      vtkSmartPointer<vtkHexahedron>        hex       =   vtkSmartPointer<vtkHexahedron>::New();

      EPOINT  EP1, EP2;

      scalars->SetNumberOfTuples(count);

      int     ngbf, ee, *tt, rr, cc, start1=0, start2=0;
      float   val;
      double  umn = 0.0, umx = 0.0;
      
      for(iii=0;iii<Npatch;iii++)
      {
         ngbf  = SolidResult[iii].ngbf;
         start2 += SolidResult[iii].nelem;
         
         SparseMatrixXd  coeffMatTmp(ngbf, ngbf);
         VectorXd  rhsVec(ngbf), soln(ngbf);
         rhsVec.setZero();

         for(ee=start1;ee<start2;ee++)
           elem[ee]->toPostprocess(vartype, varindex, index, coeffMatTmp, rhsVec);

         MatrixXd   coeffMat(coeffMatTmp);

         soln = coeffMat.ldlt().solve(rhsVec);

         ListArray<double> tmpVec;  tmpVec.setDim(ngbf);
  
         for(ii=0;ii<ngbf;ii++)
           tmpVec[ii] = soln(ii);
      
         SolidResult[iii].Values[0] = tmpVec;
         
         start1 += SolidResult[iii].nelem;
      }

//      cout << " AAAAAAAAAA " << endl;

      count=0; 
      for(iii=0;iii<Npatch;iii++)
      {
          for(kk=0;kk<n3[iii];kk++)
          {
             for(jj=0;jj<n2[iii];jj++)
             {
                 for(ii=0;ii<n1[iii];ii++)
                 {
                     EP2 = SolidResult[iii].SolidPoint(uu[iii][ii], vv[iii][jj], ww[iii][kk]).CalcEuclid();

                     points->InsertNextPoint(EP2.x, EP2.y, EP2.z);

                     val = SolidResult[iii].computeValue(1, uu[iii][ii], vv[iii][jj], ww[iii][kk]) ;              

                     if(val < umn) umn = val;
                     if(val > umx) umx = val;
          
                     scalars->SetTuple1(count++, val);
                  }
              }
          }
      }

//      cout << " AAAAAAAAAA " << endl;

      if(umnxflag)
      {
        umn = umin;
        umx = umax;
      }
 
      cout << " umin & umax ... : " << umn << '\t' << umx << endl;
      count = 0;
      for(iii=0;iii<Npatch;iii++)
      {    
          for(kk=0;kk<n3[iii]-1;kk++)
          {
              ind5 = count + n1[iii]*n2[iii]*kk;
              ind6 = count + n1[iii]*n2[iii]*(kk+1);

              for(jj=0;jj<n2[iii]-1;jj++)
              {
                 ind1 = ind5 + n1[iii]*jj;
                 ind2 = ind5 + n1[iii]*(jj+1);
       
                 ind3 = ind6 + n1[iii]*jj;
                 ind4 = ind6 + n1[iii]*(jj+1);

                 for(ii=0;ii<n1[iii]-1;ii++)
                 {
                    pts[0] = ind1+ii;          pts[4] = ind3+ii;
                    pts[1] = pts[0]+1;         pts[5] = pts[4]+1;
                    pts[3] = ind2+ii;          pts[7] = ind4+ii;
                    pts[2] = pts[3]+1;         pts[6] = pts[7]+1;
              
                    for(ll=0;ll<8;ll++)
                       hex->GetPointIds()->SetId(ll, pts[ll]);
          
                    plotvtk.uGrid->InsertNextCell(hex->GetCellType(), hex->GetPointIds());
                 }
              }
          }
          count += nn[iii];
      }
      
      
      if(Npatch > 1)
      {
          int patchnum1, patchnum2, side1, side2, nn1, nn2, iii;
          float val1, val2, valtmp;

          for(iii=0;iii<intfdata.n;iii++)//          for(iii=0;iii<1;iii++)
          {
             start1 = start2 = 0;
             patchnum1 = intfdata[iii][0] - 1;
             side1 = intfdata[iii][1];
             patchnum2 = intfdata[iii][2] - 1;
             side2 = intfdata[iii][3];
             
             for(ii=0;ii<patchnum1;ii++)
                start1 += nn[ii];
             for(ii=0;ii<patchnum2;ii++)
                start2 += nn[ii];
             
             nn1 = n1[patchnum1] * n2[patchnum1];
             nn2 = n1[patchnum2] * n2[patchnum2];
             //cout << start1 << '\t' << start2 << '\t' << nn1 << '\t' << nn2 << endl;

             if(side1 == 1)
             {
                ind1 = start1;
                ind3 = start2 + nn2*(n3[patchnum2]-1);

                for(jj=0;jj<n2[patchnum1];jj++)
                {
                   for(ii=0;ii<n1[patchnum1];ii++)
                   {
                      ind2 = ind1++;
                      ind4 = ind3++;
                   
                      //cout << ind2 << '\t' << ind4 << endl;
                   
                      val1 = scalars->GetValue(ind2);
                      val2 = scalars->GetValue(ind4);
                   
                      valtmp = 0.5*(val1 + val2);
                   
                      scalars->SetValue(ind2, valtmp);
                      scalars->SetValue(ind4, valtmp);
                   }
                }
             }

             if(side1 == 3)
             {
                for(jj=0;jj<n3[patchnum1];jj++)
                {
                   ind1 = start1 + nn1 * jj;
                   ind3 = start2 + nn2 * (jj+1) - n1[patchnum2] ;
                   for(ii=0;ii<n1[patchnum1];ii++)
                   {
                      ind2 = ind1+ii;
                      ind4 = ind3+ii;

                      //cout << ind2 << '\t' << ind4 << endl;
                   
                      val1 = scalars->GetValue(ind2);
                      val2 = scalars->GetValue(ind4);
                   
                      valtmp = 0.5*(val1 + val2);
                   
                      scalars->SetValue(ind2, valtmp);
                      scalars->SetValue(ind4, valtmp);
                   }
                }
             }
             

             if(side1 == 5)
             {
                for(jj=0;jj<n2[patchnum1];jj++)
                {
                   ind1 = start1 + n1[patchnum1]*jj;
                   ind3 = start2 + n1[patchnum2]*(jj+1)-1;
                   for(ii=0;ii<n3[patchnum1];ii++)
                   {
                      ind2 = ind1 + nn1*ii;
                      ind4 = ind3 + nn2*ii;
                   
                      //cout << ind2 << '\t' << ind4 << endl;
                   
                      val1 = scalars->GetValue(ind2);
                      val2 = scalars->GetValue(ind4);
                   
                      valtmp = 0.5*(val1 + val2);
                   
                      scalars->SetValue(ind2, valtmp);
                      scalars->SetValue(ind4, valtmp);
                   }
                }
             }
          }
      }

    mapper1->ScalarVisibilityOn();
    mapper1->SetScalarModeToUsePointData();
    mapper1->SetColorModeToMapScalars();
    mapper1->SetScalarRange(umn, umx);
    mapper1->InterpolateScalarsBeforeMappingOn();

    hueLut->SetHueRange(0.66667, 0.0);
    hueLut->SetTableRange(umn, umx);
    hueLut->SetNumberOfColors(nCol);
    hueLut->SetRampToLinear();
    hueLut->SetScaleToLinear();
    hueLut->Build();
    
    vtkSmartPointer<vtkTextProperty>    txtprop   =   vtkSmartPointer<vtkTextProperty>::New();
    
    txtprop->SetColor(0.0,0.0,0.0);
    txtprop->BoldOn();
    txtprop->ItalicOff();


    mapper1->SetLookupTable(hueLut);

    scalarBar->SetLookupTable(hueLut);
    scalarBar->SetMaximumNumberOfColors(nCol);
    scalarBar->SetNumberOfLabels(nCol+1);
    scalarBar->SetTitleTextProperty(txtprop);
    scalarBar->SetLabelTextProperty(txtprop);
    scalarBar->SetWidth(0.12);
    scalarBar->SetLabelFormat("%10.4f");
    
    plotvtk.uGrid->SetPoints(points);
    plotvtk.uGrid->GetPointData()->SetScalars(scalars);

    vtkSmartPointer<vtkXMLUnstructuredGridWriter>  writerUGridVTK =  vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    cout << " llllllllllll " << endl;    
    char fname[50];

    //sprintf(fname,"%s%06d%s", "BSpline-",filecount, ".vtu");
    sprintf(fname,"%s%s%06d%s", files.Ofile.asCharArray(),"-",filecount, ".vtu");

    writerUGridVTK->SetFileName(fname);

#if VTK_MAJOR_VERSION == 5
    writerUGridVTK->SetInput(plotvtk.uGrid);
    mapper1->SetInputConnection(plotvtk.uGrid->GetProducerPort());
#else
    writerUGridVTK->SetInputData(plotvtk.uGrid);
    mapper1->SetInputData(plotvtk.uGrid);
#endif

    writerUGridVTK->Write();

    actor1->SetMapper(mapper1);

    plotvtk.rendr->AddActor(actor1);
    plotvtk.rendr->AddActor2D(scalarBar);
    plotvtk.rendr->ResetCamera();
    plotvtk.renWindow->Render();

  return;
}






void  IsogeometricFEM::postProcessFlow(int vartype, int vardir, int nCol, bool umnxflag, double umin, double umax, int* resln)
{
   if(!plotvtk.ActiveFlag)
     plotvtk.set();

//    plotvtk.clearWindow();

    vtkSmartPointer<vtkDataSetMapper>        mapper   =  vtkSmartPointer<vtkDataSetMapper>::New();
    vtkSmartPointer<vtkActor>                actor    =  vtkSmartPointer<vtkActor>::New();
    vtkSmartPointer<vtkUnstructuredGrid>     uGrid    =  vtkSmartPointer<vtkUnstructuredGrid>::New(); 
    vtkSmartPointer<vtkPoints>               points   =  vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkQuad>                  quad    =  vtkSmartPointer<vtkQuad>::New();
    vtkSmartPointer<vtkDoubleArray>        vectors    =  vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray>        scalars    =  vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray>        vectors2   =  vtkSmartPointer<vtkDoubleArray>::New();

    int  dd, ii, jj, kk, ll, count, index, ind1, ind2, n1, n2;

    double   fact[4], geom[2], incr1, incr2, start1, start2, val, *tmp1, *tmp2;

    double   vec[3];

    VectorArray<double>  U, V, uu, vv;
    vector<vector<double> >  outp2;
    vector<double>  xx, yy;
    
    outp2.resize(ndf);

    //create_vector2(SurfaceListFinal[0].U, resln[0], uu);
    //create_vector2(SurfaceListFinal[0].V, resln[1], vv);

    create_vector(0.0,1.0,0.01,uu);
    vv = uu;
    
    EPOINT EP1, EP2;

    for(jj=0;jj<vv.n;jj++)
    {
        for(ii=0;ii<uu.n;ii++)
        {
           EP1 = SurfaceListFinal[0].SurfacePoint(uu[ii], vv[jj]).CalcEuclid();
           //EP2 = SurfaceResult[0].SurfacePoint(uu[ii], vv[jj]).CalcEuclid();

           //outp2[0].push_back(EP2.x-EP1.x);
           //outp2[1].push_back(EP2.y-EP1.y);
           
           for(kk=0;kk<ndf;kk++)
             outp2[kk].push_back(SurfaceResult[0].computeValue((kk+1),uu[ii], vv[jj]));

           xx.push_back(EP1.x);
           yy.push_back(EP1.y);

           points->InsertNextPoint(EP1.x, EP1.y, 0.0);
        }
    }
       
      // prepare and write a file to postprocess in matplotlib

      ofstream fout("post-process.dat");

      if(fout.fail())
      {
         cout << " Could not open the Output file" << endl;
      exit(1);
      }

      fout.setf(ios::fixed);
      fout.setf(ios::showpoint);
      fout.precision(8);
      
      fout << uu.n << '\t' << vv.n ;
      for(dd=0;dd<ndf;dd++)
        fout << '\t' << 0.0 ;
      fout << endl;

      for(ii=0;ii<xx.size();ii++)
      {
         fout << xx[ii] << '\t' <<  yy[ii] ;
         for(dd=0;dd<ndf;dd++)
           fout << '\t' << outp2[dd][ii] ;
         fout << endl;
      }

      fout.close();
      
      
    // create uGrid to postprocess in paraview/mayavi

    vtkIdType pt[4];
    
              // create the connectivity of elements/nodes
              // use triangle/quadrilateral as required
              
              n1 = uu.n;
              n2 = vv.n;
              
              for(jj=0;jj<vv.n-1;jj++)
              {
                 ind1 = n1 * jj ;
                 ind2 = n2 * (jj+1) ;
                 
                 for(ii=0;ii<uu.n-1;ii++)
                 {
                    pt[0] = ind1 + ii;

                    pt[1] = pt[0] + 1;

                    pt[3] = ind2 + ii;

                    pt[2] = pt[3] + 1;

                    for(ll=0;ll<4;ll++)
                       quad->GetPointIds()->SetId(ll, pt[ll]);
          
                    uGrid->InsertNextCell(quad->GetCellType(), quad->GetPointIds());
                 }
              }

    vectors->SetName("vel");
    scalars->SetName("pres");
    vectors2->SetName("vort");
    
    ii = outp2[0].size();
    
    vectors->SetNumberOfComponents(3);
    vectors->SetNumberOfTuples(ii);
    scalars->SetNumberOfTuples(ii);
        
    vec[2] = 0.0;
    
    for(ii=0;ii<outp2[0].size();ii++)
    {
      vec[0] = outp2[0][ii];
      vec[1] = outp2[1][ii];
      vectors->InsertTuple(ii, vec);
      scalars->SetTuple1(ii, outp2[2][ii]);
    }
    
    //assign nodal coordinates and field data to uGrid
    // no need to create lookup table here. All this stuff can be done in Paraview

    uGrid->SetPoints(points);

    uGrid->GetPointData()->SetScalars(scalars);
    uGrid->GetPointData()->SetVectors(vectors);
    // create a write object and write uGrid to it
    
    if(ndf == 4)
    {
       vectors2->SetNumberOfComponents(3);
       vectors2->SetNumberOfTuples(outp2[0].size());
    
       vec[0] = vec[1] = 0.0;
    
       for(ii=0;ii<outp2[0].size();ii++)
       {
          vec[2] = outp2[3][ii];
          vectors2->InsertTuple(ii, vec);
       }

       uGrid->GetPointData()->AddArray(vectors2);
    }

    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =  vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName("Bspline.vtu");

#if VTK_MAJOR_VERSION == 5
    writer->SetInput(uGrid);
    mapper->SetInputConnection(uGrid->GetProducerPort());
#else
    writer->SetInputData(uGrid);
    mapper->SetInputData(plotvtk.uGrid);
#endif

    writer->Write();

    actor->SetMapper(mapper);
    actor->GetProperty()->EdgeVisibilityOff();
    actor->GetProperty()->SetEdgeColor(0,0,0);
    actor->GetProperty()->SetPointSize(5.0);
    actor->GetProperty()->SetLineWidth(2.0);
    //actor->GetProperty()->SetColor(0,0,1);

    plotvtk.rendr->AddActor(actor);

    plotvtk.rendr->ResetCamera();
    plotvtk.renWindow->Render();


   return;
}







void  IsogeometricFEM::ConvectionSolver(int Niter, double stol, double tol, double tol2)
{
   cout << " Convection Solver ........... " << endl;
   cout << " Niter " << Niter << endl;

   int  ee, start1,  start2, iter, ii, jj, iii, kk, nn;
   start1 = start2  =  0;
   
   double  fact = 1.0, *du;

   soln.zero();
   
   firstIter = false;

   for(iter=0;iter<Niter;iter++)
   {
      solverEigen->zeroMtx();
      solverEigen->rhsVec.setZero();
      reac.zero();

      for(ee=0;ee<totnumel;ee++)
      {
         localStiffnessError = elem[ee]->calcStiffnessAndResidual();
      
         elem[ee]->AssembleElementMatrix(1, solverEigen->mtx);
         elem[ee]->AssembleElementVector(1, 0, &(solverEigen->rhsVec[0]), &(reac[0]), start1, start2);
      }

      rNorm  = solverEigen->rhsVec.norm();
      
      if(iter == 0)
      {
         for(ii=0;ii<patch.n;ii++)
           NurbsBaseResult[ii]->addInitDOFvalues();

         for(kk=0;kk<ntotgbf;kk++)
           NurbsBaseResult[0]->Values[0][kk]  = mpapTime.dt * solnInit[ii];
      }

      //cout << " rhsVec " << endl;        printVector(&(rhsVec[0]), rhsVec.n);

      COUT << domain.name(this); printf("  %11.4e\n",rNorm);

      solverEigen->currentStatus = ASSEMBLY_OK;

      solverEigen->factoriseAndSolve();

      //cout << " result " << endl;        printVector(du, ntoteqs);

      for(kk=0;kk<ntoteqs1;kk++)
        soln[assy4r[kk]] = solverEigen->soln[kk];

      for(kk=0;kk<ntotgbf;kk++)
        NurbsBaseResult[0]->Values[0][kk]  = soln[kk];
  }

  for(iii=0;iii<patch.n;iii++)
    //NurbsBaseResult[iii]->updateCoordinates(&(soln[0]));
    NurbsBaseResult[iii]->updateCoordinates(&(NurbsBaseResult[iii]->Values[0][0]));

  return;
}






void  IsogeometricFEM::LMSolver(int Niter, double stol, double tol, double tol2)
{
   cout << " LM Solver ........... " << endl;
   cout << " Niter " << Niter << endl;

   int  e, start1,  start2, iter, ii, jj, iii, kk, nn, temp;
   start1 = start2  =  0;

   double  fact = 0.001, *du;

   soln.zero();
   
   firstIter = false;

   for(iii=0;iii<Npatch;iii++)
   NurbsBaseResult[iii]->geomToVector(&(solnPrev[0]));

   for(ii=0;ii<ndf;ii++)
     Values[ii] = NurbsBaseResult[0]->Values[ii] ;

   for(iter=0;iter<Niter;iter++)
   {
      solverEigen->zeroMtx();
      solverEigen->rhsVec.setZero();
      reac.zero();

      for(e=0;e<totnumel;e++)
      {
         localStiffnessError = elem[e]->calcStiffnessAndResidual();
      
         elem[e]->AssembleElementMatrix(1, solverEigen->mtx);
         elem[e]->AssembleElementVector(firstIter, 0, &(solverEigen->rhsVec[0]), &(reac[0]), start1, start2);
      }

      for(e=0;e<totnumel;e++)  // loop over all the elements
      {
        if(elem[e]->tracflag)
        {
           temp = elem[e]->calcLoadVector();
           elem[e]->AssembleElementMatrix(1, solverEigen->mtx);
           elem[e]->AssembleElementVector(firstIter, 0, &(solverEigen->rhsVec[0]), &(reac[0]), 0, 0);
        }
      }

      //cout << " rhsVec " << endl;        printVector(&(rhsVec[0]), rhsVec.n);

      rNormPrev = rNorm;
      rNorm     = solverEigen->rhsVec.norm();
      
      COUT << domain.name(this); printf("  %12.6f \t %11.4e \t %11.4e\n", fact, rNormPrev, rNorm);

      if(converged())
        break;

      for(e=0;e<totnumel;e++)
        elem[e]->AssembleElementMatrix3(1, fact, solverEigen->mtx);

      solverEigen->currentStatus = ASSEMBLY_OK;

      solverEigen->factoriseAndSolve();

      for(kk=0;kk<ntoteqs1;kk++)
        soln[assy4r[kk]]  = solverEigen->soln[kk];

      for(kk=0;kk<ntotgbf;kk++)
      {
         iii = kk*ndf;
         for(ii=0;ii<ndf;ii++)
           NurbsBaseResult[0]->Values[ii][kk]  += soln[iii+ii];
      }

      for(iii=0;iii<Npatch;iii++)
        NurbsBaseResult[iii]->updateCoordinates(&(soln[0]));
      
      calcAndAssyInternalForceVector(1.0);

      rNormPrev = rNorm;
      rNorm     = solverEigen->rhsVec.norm();
      
      printf("\t%12.6e\t%12.6e\n",rNormPrev, rNorm);

      while(rNorm > rNormPrev)
      {
         for(iii=0;iii<patch.n;iii++)
           NurbsBaseResult[iii]->resetGeometry(&(solnPrev[0]));

         for(kk=0;kk<ndf;kk++)
           NurbsBaseResult[0]->Values[kk] = Values[kk];

         //for(e=0;e<totnumel;e++)
           //elem[e]->AssembleElementMatrix3(1, -fact, solverEigen->mtx);

         fact *= 10.0;

         for(e=0;e<totnumel;e++)
           elem[e]->AssembleElementMatrix3(1, fact, solverEigen->mtx);

         //calcAndAssyInternalForceVector(1.0);

         solverEigen->currentStatus = ASSEMBLY_OK;

         solverEigen->factoriseAndSolve();

         for(kk=0;kk<ntoteqs1;kk++)
           soln[assy4r[kk]]  = solverEigen->soln[kk];

         for(kk=0;kk<ntotgbf;kk++)
         {
            iii = kk*ndf;
            for(ii=0;ii<ndf;ii++)
              NurbsBaseResult[0]->Values[ii][kk]  += soln[iii+ii];
         }

        for(iii=0;iii<Npatch;iii++)
           NurbsBaseResult[iii]->updateCoordinates(&(soln[0]));

         calcAndAssyInternalForceVector(1.0);

         rNorm  = solverEigen->rhsVec.norm();

         //printf(" \t \t \t %11.4e\n",rNorm);
         printf("  %12.6f \t %11.4e\n", fact, rNorm);
      }

      fact /= 50.0;

      for(iii=0;iii<Npatch;iii++)
        NurbsBaseResult[iii]->geomToVector(&(solnPrev[0]));

      for(ii=0;ii<ndf;ii++)
        Values[ii] = NurbsBaseResult[0]->Values[ii] ;
  }

  return;
}




