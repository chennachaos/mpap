#include <Eigen/Dense>

#include <math.h>
#include "Debug.h"
#include "MpapTime.h"
#include "PlotVTK.h"
#include "NurbsElem3DStructSolid.h"
#include "NurbsUtilitiesSOLID.h"
#include <assert.h>
#include "ComputerTime.h"


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

using namespace std;
using namespace Eigen;


extern ComputerTime       computerTime;
extern MpapTime mpapTime;
extern PlotVTK plotvtk;




NurbsElem3DStructSolid::NurbsElem3DStructSolid(void)
{
  if (debug) cout << " constructor NurbsElem3DStructSolid\n\n";
}



NurbsElem3DStructSolid::~NurbsElem3DStructSolid()
{
  if (debug) cout << " destructor NurbsElem3DStructSolid\n\n";
}

/*
       cout << '\t' << " deformation Gradient  " << endl;
       ttt=0;
        for(ii=0;ii<3;ii++)
        {
           printf("\t");
           for(jj=0;jj<3;jj++)
             printf("%14.10f\t", F[ttt++]*100000);
           printf("\n");
        }
        cout << endl;
        cout << endl;
*/


/*
int NurbsElem3DStructSolid::calcStiffnessAndResidual()
{
   char fct[] = "NurbsElem3DStructSolid::calcStiffnessAndResidual";

   computerTime.go(fct);
   
   MatrixXd  F(3,3), cc(6,6), bc(3,6), stifftmp(nsize,nsize), Bmat(6, nsize);
   bc.setZero();
   stifftmp.setZero();
   Bmat.setZero();
   
   VectorXd  stre(6);

   double  detF=0.0, fact, fact1, fact2, fact3, fact4;

   double  dvol, dvol0, Jac, dt, totvol, bb1, bb2, bb3, dN_dx[nlbf], dN_dy[nlbf], dN_dz[nlbf];

   int  index, ll = 0, ii, jj, gp1, gp2, gp3, kk, threeI, threeIp1, threeIp2, threeJ, threeJp1, threeJp2;

   int  err = 0,  isw = 3,  count = 1, count1 = 0, ttt;

   for(ii=0;ii<nsize;ii++)
      stiffness_local[ii].zero();
   resi.zero();

   double *gaussweights = &(solid0->gaussweights[0]);

   totvol = 0.0;
   dt = mpapTime.dt;

     count1 = 0;
     for(gp3=0;gp3<nGP3;gp3++)
     {
     for(gp2=0;gp2<nGP2;gp2++)
     {
     for(gp1=0;gp1<nGP1;gp1++)
     {
             index = count1*3;
             
             solid0->ShapeFunDerivatives(&startindex[0], &(knotsAtGPs[index]), dN_dx, dN_dy, dN_dz, Jac);
             
             fact = gaussweights[count1] * JacMultFact;
             dvol0 = Jac * fact;
             dvol  = dvol0;
             
             //totvol += dvol0;

             solid1->deformationGradient(&startindex[0], 1, dN_dx, dN_dy, dN_dz, &F(0,0), detF);

             if(detF < 0.0)   return 1;

             if(finite) 
             {
                solid1->ShapeFunDerivatives(&startindex[0], &(knotsAtGPs[index]), dN_dx, dN_dy, dN_dz, Jac);
                dvol = Jac * fact;
             }             
             //printf("Jac... %12.8f\n", Jac);

             matlib3d_(matDat, &F(0,0), &stre(0), &cc(0,0), &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &isw, &err, &count, NULL);
             count++;
             count1++;
             ll += nivGP;

             if(err !=0)
               return 1;
               
              stre *= dvol;
              cc   *= dvol;


             //printf(" stresses \n");        printf("\t%12.8f\t%12.8f\t%12.8f\t%12.8f\t%12.8f\t%12.8f\n\n", stre[0], stre[1], stre[2], stre[3], stre[4], stre[5]);

             //==============================================
             // CALCULATE TANGENT STIFFNESS
             //==============================================

            //   part 1. -- material part (not necessarily symmetric!!)

            for(ii=0;ii<nlbf;ii++)
            {
                threeI   = 3*ii;
                threeIp1 = threeI+1;
                threeIp2 = threeI+2;

                Bmat(0,threeI)    =  dN_dx[ii];
                Bmat(1,threeIp1)  =  dN_dy[ii];
                Bmat(2,threeIp2)  =  dN_dz[ii];
                
                Bmat(3,threeI)   =  dN_dy[ii]; Bmat(3,threeIp1)  =  dN_dx[ii];
                Bmat(4,threeIp1) =  dN_dz[ii]; Bmat(4,threeIp2)  =  dN_dy[ii];
                Bmat(5,threeI)   =  dN_dz[ii]; Bmat(5,threeIp2)  =  dN_dx[ii];
                
            }

            stifftmp += Bmat.transpose() * cc * Bmat;

     } //gp1
     } //gp2
     } //gp3

//  cout << '\t' << " Total Volume for element # " << elenum << " is = " << totvol << endl; cout << endl;

   computerTime.stopAndPrint(fct);

 // if(elenum == 0)     printMatrix(stifftmp);

//cout << patchnum << '\t' << elenum << endl;
//printForceVector();

  return 0;
}
*/



int NurbsElem3DStructSolid::calcStiffnessAndResidual()
{

  double  F[9], detF=0.0, fact, fact1, fact2, fact3, stre[6], cc[6][6], bc[3][6], bbf[3];
  double  dvol, dvol0, Jac, dt, totvol, bb1, bb2, bb3, N[nlbf], dN_dx[nlbf], dN_dy[nlbf], dN_dz[nlbf];

  int  index, ll = 0, ii, jj, gp1, gp2, gp3, kk, TI, TIp1, TIp2, TJ, TJp1, TJp2;
  int  err = 0,  isw = 3,  count = 1, count1 = 0, ttt;

  Klocal.setZero();
  Flocal.setZero();

  double *gaussweights = &(solid0->gaussweights[0]);

    totvol = 0.0;
    dt = mpapTime.dt;

    count1 = 0;
    for(gp3=0;gp3<nGP3;gp3++)
    {
    for(gp2=0;gp2<nGP2;gp2++)
    {
    for(gp1=0;gp1<nGP1;gp1++)
    {
        index = count1*3;

        //solid0->ShapeFunctions(knotsAtGPs[index], knotsAtGPs[index+1], knotsAtGPs[index+2], NN);

        solid0->ShapeFunDerivatives(&startindex[0], &(knotsAtGPs[index]), N, dN_dx, dN_dy, dN_dz, Jac);

        fact = gaussweights[count1] * JacMultFact;
        dvol0 = Jac * fact;
        dvol  = dvol0;

        solid1->deformationGradient(&startindex[0], 1, dN_dx, dN_dy, dN_dz, F, detF);

        if(detF < 0.0)   return 1;

        if(finite) 
        {
          solid1->ShapeFunDerivatives(&startindex[0], &(knotsAtGPs[index]), N, dN_dx, dN_dy, dN_dz, Jac);
          dvol = Jac * fact;
        }             
        //printf("Jac... %12.8f\n", Jac);

        matlib3d_(matDat, F, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &isw, &err, &count, NULL);
        count++;
        count1++;
        ll += nivGP;

        if(err !=0)
          return 1;

        for(ii=0;ii<6;ii++)
        {
          stre[ii] *= dvol;
          for(jj=0;jj<6;jj++)
            cc[ii][jj] *= dvol;
        }

        //printf(" stresses \n");
        //printf("\t%20.18f\t%20.18f\t%20.18f\t%20.18f\t%20.18f", F[0], F[1], F[2], F[3], F[4]);
        //printf("\t%20.18f\t%20.18f\t%20.18f\t%20.18f\n\n", F[5], F[6], F[7], F[8]);
        //printf("\t%14.12f\t%14.12f\t%14.12f\t%14.12f\t%14.12f\t%14.12f\n\n", stre[0], stre[1], stre[2], stre[3], stre[4], stre[5]);

        fact  = rho0 * dvol0 ;
        bbf[0] = bforce[0] * fact ;
        bbf[1] = bforce[1] * fact ;
        bbf[2] = bforce[2] * fact ;
        bbf[0] = bbf[1] = bbf[2] = 0.0;

        //==============================================
        // CALCULATE TANGENT STIFFNESS
        //==============================================

        //   part 1. -- material part (not necessarily symmetric!!)

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

          Flocal[TI]   += (N[ii] * bbf[0] - bb1*stre[0] - bb2*stre[3] - bb3*stre[5]) ;
          Flocal[TIp1] += (N[ii] * bbf[1] - bb1*stre[3] - bb2*stre[1] - bb3*stre[4]) ;
          Flocal[TIp2] += (N[ii] * bbf[2] - bb1*stre[5] - bb2*stre[4] - bb3*stre[2]) ;

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
        }

        //   part 2. -- geometrical matrix  (if geometry nonlinear)

        if(finite)
        {
          for(ii=0;ii<nlbf;ii++)
          {
            bb1 = dN_dx[ii];
            bb2 = dN_dy[ii];
            bb3 = dN_dz[ii];

            fact1 = bb1*stre[0] + bb2*stre[3] + bb3*stre[5] ;
            fact2 = bb1*stre[3] + bb2*stre[1] + bb3*stre[4] ;
            fact3 = bb1*stre[5] + bb2*stre[4] + bb3*stre[2] ;

            TJ = 3*ii;
            for(jj=0;jj<nlbf;jj++)
            {
              fact = fact1 * dN_dx[jj] + fact2 * dN_dy[jj] + fact3 * dN_dz[jj];

              TJ = 3*jj;

              Klocal(TI,   TJ)   += fact ;
              Klocal(TI+1, TJ+1) += fact ;
              Klocal(TI+2, TJ+2) += fact ;
            }
          }
        } //if(finite)
    } //gp1
    } //gp2
    } //gp3
    
    //printForceVector();

  return 0;
}




int NurbsElem3DStructSolid::calcInternalForces()
{
/*
//   computerTime.go(fct);

   double  F[9], detF=0.0, fact, fact1, fact2, fact3, fact4, stre[6], cc[6][6];

   double  dvol, dvol0, Jac, dt, totvol, bb1, bb2, bb3, N[nlbf], dN_dx[nlbf], dN_dy[nlbf], dN_dz[nlbf];

   int  index, ll = 0, ii, jj, gp1, gp2, gp3, kk, threeI, threeIp1, threeIp2;

   int  err = 0,  isw = 3,  count = 1, count1 = 0, ttt;

   resi.zero();

   double *gaussweights = &(solid0->gaussweights[0]);

   totvol = 0.0;
   dt = mpapTime.dt;

     count1 = 0;
     for(gp3=0;gp3<nGP3;gp3++)
     {
     for(gp2=0;gp2<nGP2;gp2++)
     {
     for(gp1=0;gp1<nGP1;gp1++)
     {
             index = count1*3;
             
             solid0->ShapeFunDerivatives(&startindex[0], &(knotsAtGPs[index]), N, dN_dx, dN_dy, dN_dz, Jac);
             
             fact = gaussweights[count1] * JacMultFact;
             dvol0 = Jac * fact;
             dvol  = dvol0;
             
             //totvol += dvol0;

             solid1->deformationGradient(&startindex[0], 1, dN_dx, dN_dy, dN_dz, F, detF);

             if(detF < 0.0)   return 1;

             if(finite) 
             {
                solid1->ShapeFunDerivatives(&startindex[0], &(knotsAtGPs[index]), N, dN_dx, dN_dy, dN_dz, Jac);
                dvol = Jac * fact;
             }             
             //printf("Jac... %12.8f\n", Jac);

             matlib3d_(matDat, F, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &isw, &err, &count, NULL);
             count++;
             count1++;
             ll += nivGP;

             for(ii=0;ii<6;ii++)
               stre[ii] *= dvol;

             //printf(" stresses \n");        printf("\t%12.8f\t%12.8f\t%12.8f\t%12.8f\t%12.8f\t%12.8f\n\n", stre[0], stre[1], stre[2], stre[3], stre[4], stre[5]);

             //==============================================
             // CALCULATE TANGENT STIFFNESS
             //==============================================

            //   part 1. -- material part (not necessarily symmetric!!)

            for(ii=0;ii<nlbf;ii++)
            {
                bb1 = dN_dx[ii];
                bb2 = dN_dy[ii];
                bb3 = dN_dz[ii];
                
                threeI   = 3*ii;
                threeIp1 = threeI+1;
                threeIp2 = threeI+2;

                resi[threeI]   -= (bb1*stre[0] + bb2*stre[3] + bb3*stre[5]) ;
                resi[threeIp1] -= (bb1*stre[3] + bb2*stre[1] + bb3*stre[4]) ;
                resi[threeIp2] -= (bb1*stre[5] + bb2*stre[4] + bb3*stre[2]) ;
                
             }

     } //gp1
     } //gp2
     } //gp3
  
*/
//printForceVector();

  return 0;
}



void NurbsElem3DStructSolid::discreteContourplot(int vartype, int varindex, int index, int nCol, double umin, double umax)
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


    vtkIdType pts[8], pt1, pt2;

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




void NurbsElem3DStructSolid::projectToKnots(bool extrapolateFlag, int vartype, int varindex, int index)
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




void NurbsElem3DStructSolid::projectStress(int varindex, double* outval)
{
/*
    if(varindex > 8)
    {
       cout << '\t' << "    NurbsElem3DStructSolid::projectStress .... : Error in 'varindex' " << endl;
       return;
    }


   double F[9], detF, stre[6], cc[6][6], Jac, dt;

   int   err    = 0,
         isw    = 3,
         count  = 1, 
         count1 = 0, index, ll = 0, ii, jj, gp1, gp2, gp3;

   double  N[nlbf], dN_dx[nlbf], dN_dy[nlbf], dN_dz[nlbf];

   dt = mpapTime.dt;

   int nivEL = nGP * nivGP;
   for(ii=0;ii<nivEL;ii++)
     intVar2[ii] = intVar1[ii];

    // loop over Gauss points
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

          if(varindex < 6)
             outval[count1] = stre[varindex];
          else if(varindex == 6)
             outval[count1] = vonMises3D(stre);
          else if(varindex == 7)
             outval[count1] = (stre[0]+stre[1]+stre[2])/3.0;

          count++;
          count1++;
          ll += nivGP;

    }//gp1
    }//gp2
    }//gp3
*/
   return;
}



void NurbsElem3DStructSolid::projectStrain(int vartype, int varindex, double* outval)
{
    cout << "        NurbsElem3DStructSolid::projectStrain.......Not implemented yet   " << endl;
 
  return;
}




void NurbsElem3DStructSolid::projectIntVar(int index, double* outval)
{
   int ind1, ii, jj, kk;

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






int NurbsElem3DStructSolid::calcMassMatrix(int lumpInd, double dt)
{
    cout << "        NurbsElem3DStructSolid::calcMassMatrix.......Not implemented yet   " << endl;
  return 0;
}





int NurbsElem3DStructSolid::calcOutput(double u1, double v1)
{
    cout << "        NurbsElem3DStructSolid::calcOutput.......Not implemented yet   " << endl;

  return 0;
}




void NurbsElem3DStructSolid::toPostprocess(int vartype, int varindex, int type, SparseMatrixXd&  coeffMat, VectorXd& rhsVec)
{
/*
//      cout << " SSSSSSSSSSS " << elenum << endl;

   MatrixXd  Nlocal(nlbf,nlbf);
   VectorXd  NN(nlbf), rhslocal(nlbf);

   Nlocal.setZero();   
   rhslocal.setZero();


   double F[9], detF, stre[6], cc[6][6], Jac, dt, pres, utemp, vtemp, wtemp;

   int   err    = 0,
         isw    = 3,
         count  = 1, 
         count1 = 0, index, ll = 0, gp1, gp2, gp3, ii, jj, row, col;

   double  N[nlbf], dN_dx[nlbf], dN_dy[nlbf], dN_dz[nlbf];

   dt = mpapTime.dt;

   int nivEL = nGP * nivGP;
   for(ll=0;ll<nivEL;ll++)
     intVar2[ll] = intVar1[ll];


//      cout << " SSSSSSSSSSSSSS " << elenum << endl;


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

          Nlocal += NN * NN.transpose();
        
          if(varindex < 6)
             rhslocal += NN * stre[varindex];
          else if(varindex == 6)
             rhslocal += NN * vonMises3D(stre);
          else if(varindex == 7)
          {
             pres = solid2->computeValue(1, utemp, vtemp, wtemp);
             rhslocal += NN * pres;
          }

          count++;
          count1++;
          ll += nivGP;

    }//gp1
    }//gp2
    }//gp3


//      cout << " DDDDDDDDDDD " << elenum << endl;

    int *tt;
    tt = &(solid0->IEN[elenum][0]);
      
    for(ii=0;ii<nlbf;ii++)
    {
       row = tt[ii];
       rhsVec(row) += rhslocal(ii);
       for(jj=0;jj<nlbf;jj++)
       {
          col = tt[jj];
          coeffMat.coeffRef(row, col) += Nlocal(ii, jj);
       }
    }

//      cout << " DDDDDDDDDDD " << elenum << endl;
*/
  return;
}

/*
void  NurbsElem3DStructSolid::AssembleElementMatrix(int index, Mat mtx, int start1, int start2)
{
    PetscErrorCode ierr;
    int  ii, jj, nn=0, aa, bb, size2, ind, *tt1, *tt2;

    tt1 = &(solid0->LM[elenum][0]);

    //cout << elenum << '\t' << elenum2 << endl;
    //cout << nsize << '\t' << size2 << endl;
    //cout << surf0->LM[elenum] << endl;
    //cout << surf2->LM[elenum2] << endl;

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
      }
    }

  return;
}
*/



void  NurbsElem3DStructSolid::AssembleElementMatrix(int index, SparseMatrixXd& mtx, int start1, int start2)
{
    int  ii, jj, nn=0, aa, bb, ind, *tt1, *tt2;

    tt1 = &(solid0->LM[elenum][0]);

    /*
   printMatrix(Klocal);
   printf("\n\n\n");
   */

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
        //cout << " ii = " << ii << endl;
      }
    }


  return;
}






