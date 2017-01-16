
#include <math.h>
#include "Debug.h"
#include "MpapTime.h"
#include "Plot.h"
#include "NurbsElemMindlinPlate.h"
#include "NurbsShapeFunctions.h"
#include <assert.h>
#include "ComputerTime.h"

using namespace std;

extern ComputerTime       computerTime;
extern MpapTime mpapTime;
extern Plot plot;


NurbsElemMindlinPlate::NurbsElemMindlinPlate()
{
  if (debug) cout << " constructor NurbsElemMindlinPlate\n\n";
}



NurbsElemMindlinPlate::~NurbsElemMindlinPlate()
{
  if (debug) cout << " destructor NurbsElemMindlinPlate\n\n";

}



int NurbsElemMindlinPlate::calcStiffnessAndResidual()
{
/*
//   char fct[] = "NurbsElemMindlinPlate::calcStiffnessAndResidual";

//   computerTime.go(fct);

   double   F[4], detF=0.0, F33=0.0, fact, fact1, fact2, fact3, fact4;
   double  stre[4], cc[4][4], bc[2][4], AA, bb1, bb2, bb3, bb4;

   double   Iz, darea, temp, Jac, dt;

   int   err    = 0,
         isw    = 3,
         count  = 1,
         count1 = 0, index, ll = 0, ii, jj, gp1, gp2, row, col;

   thick = elmDat[4];

   Iz =pow(thick,3.0)/12.0;

   for(ii=0;ii<nsize;ii++)
       stiffness_local[ii].zero();

   resi.zero();

   double *gaussweights = &(surf0->gaussweights[0]);

   double   N[nlbf], dN_dx[nlbf], dN_dy[nlbf];

   // loop over Gauss points
   for(gp2=0;gp2<nGP2;gp2++)
   {
      for(gp1=0;gp1<nGP1;gp1++)
      {
         index = count1*2;

         surf0->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, Jac);

         if(Jac < 0.0)   return 1;

         darea = Jac * gaussweights[count1] ;

         F[0] = F[1] = F[2] = F[3] = 0.0;

         dt = mpapTime.dt;
         matlib2d_(matDat, F, &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, (Element*) this);
         count++;
         count1++;
         ll += nivGP;

         AA = cc[3][3] * thick * darea;

         fact = Iz * darea;

         for(ii=0;ii<4;ii++)
         {
           stre[ii] *= darea;
           for(jj=0;jj<4;jj++)
             cc[ii][jj] *= fact;
         }


         //==============================================
         // CALCULATE TANGENT STIFFNESS
         //==============================================

         for(ii=0;ii<nlbf;ii++)
         {
             bb1 = dN_dx[ii];
             bb2 = dN_dy[ii];

             bc[0][0] = (bb1 * cc[0][0] + bb2 * cc[3][0]);
             bc[0][1] = (bb1 * cc[0][1] + bb2 * cc[3][1]);
             bc[0][2] = (bb1 * cc[0][3] + bb2 * cc[3][3]);

             bc[1][0] = (bb2 * cc[1][0] + bb1 * cc[3][0]);
             bc[1][1] = (bb2 * cc[1][1] + bb1 * cc[3][1]);
             bc[1][2] = (bb2 * cc[1][3] + bb1 * cc[3][3]);

             row = ndof * ii;

             fact1 = AA * bb1;
             fact2 = AA * bb2;

             fact = AA * N[ii];

             for(jj=0;jj<nlbf;jj++)
             {
                col = ndof * jj;

                bb3 = dN_dx[jj];
                bb4 = dN_dy[jj];

                stiffness_local[row][col]     +=  AA * (bb1 * bb3 + bb2 * bb4) ;
                stiffness_local[row][col+1]   +=  fact1 * N[jj];
                stiffness_local[row][col+2]   +=  fact2 * N[jj];
                stiffness_local[row+1][col]   +=  fact * bb3;
                stiffness_local[row+2][col]   +=  fact * bb4;

                fact3 = fact * N[jj];

                stiffness_local[row+1][col+1]  +=  (bc[0][0] * bb3 + bc[0][2] * bb4 + fact3) ;
                stiffness_local[row+1][col+2]  +=  (bc[0][1] * bb4 + bc[0][2] * bb3) ;
                stiffness_local[row+2][col+1]  +=  (bc[1][0] * bb3 + bc[1][2] * bb4) ;
                stiffness_local[row+2][col+2]  +=  (bc[1][1] * bb4 + bc[1][2] * bb3 + fact3) ;
             }
         }
     //
         // internal forces
         for(ii=0;ii<nlbf;ii++)
         {
            twoI = 2*ii;
            resi[twoI]   -= (dN_dx[ii]*stre[0] + dN_dy[ii]*stre[3]) ;
            resi[twoI+1] -= (dN_dx[ii]*stre[3] + dN_dy[ii]*stre[1]) ;
         }
     //
     }
  }

//  cout << '\t' << " Total Volume for element # " << elenum << " is = " << totvol << endl; cout << endl;

//   computerTime.stopAndPrint(fct);

//printStiffnessMatrix();
//printForceVector();
*/
  return 0;
}




int NurbsElemMindlinPlate::calcMassMatrix(int lumpInd, double dt)
{
  return 0;
}




int NurbsElemMindlinPlate::calcLoadVector()
{
/*
   double  fact, Jac,  N[nlbf], dN_dx[nlbf], dN_dy[nlbf], load;

   int  count1 = 0, index, ii, jj, gp1, gp2;

   load = elmDat[5];

   resi.zero();

   double *gaussweights = &(surf0->gaussweights[0]);

   for(gp2=0;gp2<nGP2;gp2++)
   {
      for(gp1=0;gp1<nGP1;gp1++)
      {
         index = count1*2;
         count1++;

         NurbsShapeFunctions2DAlg11(surf0, startindex[0], startindex[1], knotsAtGPs[index], knotsAtGPs[index+1], N, dN_dx, dN_dy, Jac);

         if(Jac < 0.0)   return 1;

         fact = Jac * gaussweights[count1] * load;

         for(ii=0;ii<nlbf;ii++)
            resi[ndof*ii] += N[ii] * fact;
     }
  }
*/
  return 0;
}





void NurbsElemMindlinPlate::discreteContourplot(int vartype, int varindex, int index, int nCol, double umin, double umax)
{
  if(index > nivGP)
  {
     cout << '\t' << " Error in NurbsElemMindlinPlate::contourplot " << endl;
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

/*
  int tt=0;
  for(int jj=0;jj<nGP2;jj++)
  {
     for(int ii=0;ii<nGP1;ii++)
     {
        cout << '\t' << outval[tt];
        tt++;
     }
  }

  cout << endl;
  cout << endl;
*/
  double uu, vv, du, dv;

  du = uvalues[2]/nGP1;
  dv = vvalues[2]/nGP2;

  ListArray<EPOINT> S1;
  S1.setDim( (nGP1+1)*(nGP2+1) );

  int count=0, ii, jj;
  vv = vvalues[0];

  if(finite)
  {
     for(jj=0;jj<=nGP2;jj++)
     {
        uu = uvalues[0];
        for(ii=0;ii<=nGP1;ii++)
        {
           S1[count] = surf1->SurfacePoint(uu, vv).CalcEuclid();
           count++;
           uu += du;
        }
        vv += dv;
     }
  }
  else
  {
     for(jj=0;jj<=nGP2;jj++)
     {
        uu = uvalues[0];
        for(ii=0;ii<=nGP1;ii++)
        {
           S1[count] = surf0->SurfacePoint(uu, vv).CalcEuclid();
           count++;
           uu += du;
        }
        vv += dv;
     }

  }

  EPOINT *EP;

  double x1[2], x2[2], x3[2], x4[2], u1;

  int ind1, ind2, nGP1p1;

  nGP1p1 = nGP1 + 1;

  count=0;
  for(jj=0;jj<nGP2;jj++)
  {
      ind1 = nGP1p1*jj;
      ind2 = nGP1p1*(jj+1);

      for(ii=0;ii<nGP1;ii++)
      {
          u1 = outval[count];

          EP = &(S1[ind1+ii]);
          x1[0] = EP->x; x1[1] = EP->y;

          EP = &(S1[ind2+ii]);
          x2[0] = EP->x; x2[1] = EP->y;

          EP = &(S1[ind2+ii+1]);
          x3[0] = EP->x; x3[1] = EP->y;

          EP = &(S1[ind1+ii+1]);
          x4[0] = EP->x; x4[1] = EP->y;

          // contour plot for 1st triangle
          plot.triangleContourPlot(x1, x2, x3, u1, u1, u1, umin, umax, nCol);

          // contour plot for 2nd triangle
          plot.triangleContourPlot(x1, x3, x4, u1, u1, u1, umin, umax, nCol);

          count++;
      }
  }

  return;
}











void NurbsElemMindlinPlate::projectToKnots(bool extrapolateFlag, int vartype, int varindex, int index)
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


   assert(vals2project.n == 4);

    if(extrapolateFlag)
    {
       for(int ii=0;ii<4;ii++)
         vals2project[ii] = extrapolate(nGP1, (ii+1), outval);
    }
    else
    {
       vals2project[0] = outval[0];
       vals2project[1] = outval[nGP1-1];
       vals2project[2] = outval[nGP1*(nGP2-1)];
       vals2project[3] = outval[nGP1*nGP2-1];
    }

//   cout << '\t' << vals2project << endl; cout << endl;

  return;
}




void NurbsElemMindlinPlate::projectStress(int varindex, double* outval)
{

  return;
}




void NurbsElemMindlinPlate::projectStrain(int vartype, int varindex, double* outval)
{

  return;
}






void NurbsElemMindlinPlate::projectIntVar(int index, double* outval)
{
   int ind1, ii, jj;

   ind1 = 0;
   for(jj=0;jj<nGP2;jj++)
   {
       for(ii=0;ii<nGP1;ii++)
       {
           outval[ind1] = intVar2[ind1*nivGP+index];
           ind1++;
       }
   }


   return;
}



int NurbsElemMindlinPlate::calcStiffnessMatrix(double dt)
{

  return 0;
}


int NurbsElemMindlinPlate::calcInternalForces()
{

  return 0;
}



int NurbsElemMindlinPlate::calcOutput(double u1, double v1)
{

  return 0;
}




