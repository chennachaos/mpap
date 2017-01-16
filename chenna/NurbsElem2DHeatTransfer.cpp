
#include <math.h>
#include "Debug.h"
#include "MpapTime.h"
#include "Plot.h"
#include "NurbsElem2DHeatTransfer.h"
#include "NurbsShapeFunctions.h"
#include <assert.h>
#include "ComputerTime.h"
#include "Functions.h"

#include "ExactSolutionsElasticity.h"

using namespace std;
using namespace Eigen;

extern ComputerTime       computerTime;
extern MpapTime mpapTime;
extern Plot plot;


NurbsElem2DHeatTransfer::NurbsElem2DHeatTransfer()
{
  if (debug) cout << " constructor NurbsElem2DHeatTransfer\n\n";
}



NurbsElem2DHeatTransfer::~NurbsElem2DHeatTransfer()
{
  if (debug) cout << " destructor NurbsElem2DHeatTransfer\n\n";
}



/*
int NurbsElem2DHeatTransfer::calcStiffnessAndResidual()
{
   double  fact, dvol0, Jac, kx, ky, mu, s, res;

   double  N[nlbf], dN_dx[nlbf], dN_dy[nlbf], d2N_dx2[nlbf], d2N_dy2[nlbf], temp[nlbf];

   int   count1, ii, jj, gp1, gp2, index;

   kx = elmDat[4];
   ky = elmDat[5];
   mu = elmDat[6];
   
   for(ii=0;ii<nsize;ii++)
     stiffness_local[ii].zero();

   resi.zero();

   double *gaussweights = &(surf0->gaussweights[0]);
   double  *values1 = &(surf1->Values[0][0]);

   int *tt = &(surf0->IEN[elenum][0]);

   count1 = 0;
   for(gp2=0;gp2<nGP2;gp2++)
   {
   for(gp1=0;gp1<nGP1;gp1++)
   {
        index = count1*2;

        surf0->ShapeFunDerivatives2(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, d2N_dx2, d2N_dy2, Jac);

        dvol0 = Jac * gaussweights[count1] * JacMultFact;
        count1++;

        res = 0.0;
        for(ii=0;ii<nlbf;ii++)
        {
          //printf("\t %5d \t\t %12.8f \t\t %12.8f \t\t %12.8f \t\t %12.8f \n ", ii, dN_dx[ii], dN_dy[ii], d2N_dx2[ii], d2N_dy2[ii]);
          temp[ii] = d2N_dx2[ii] + d2N_dy2[ii];

          fact = values1[tt[ii]];

          res -= fact * temp[ii];
        }

        for(ii=0;ii<nlbf;ii++)
        {
           fact = temp[ii] * dvol0;

           resi[ii] += (res * fact) ;

           for(jj=0;jj<nlbf;jj++)
             stiffness_local[ii][jj]  +=  (fact * temp[jj]) ;
        }
  }//gp1
  }//gp2
  //printf("\n\n");

  return 0;
}
*/



int NurbsElem2DHeatTransfer::calcStiffnessAndResidual()
{
/*
  double F[2], dvol0, Jac, bb1, bb2, bb3, force, fact, xx, yy;
  double  N[nlbf], dN_dx[nlbf], dN_dy[nlbf];
  
  EPOINT  EP;

  int   count,  count1, index, ii, jj, gp1, gp2;
  
  LaplaceNew  analy;

  Klocal.setZero();
  Flocal.setZero();

  double *gaussweights1 = &(surf0->gaussweights1[0]);
  double *gaussweights2 = &(surf0->gaussweights2[0]);

  double  *values1 = &(surf1->Values[0][0]);

  int *tt = &(surf0->IEN[elenum][0]);

  count1 = 0;
  for(gp2=0;gp2<nGP2;gp2++)
  {
  for(gp1=0;gp1<nGP1;gp1++)
  {
      index = count1*2;

        surf0->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, Jac);
      
        fact = gaussweights2[gp2] * gaussweights1[gp1] * thick * JacMultFact;

        dvol0 = Jac * fact;
        count1++;

          EP = surf0->SurfacePoint(knotsAtGPs[index], knotsAtGPs[index+1]).CalcEuclid();
          
          force = analy.force(EP.x, EP.y);


          F[0] = F[1] = 0.0;
          for(ii=0;ii<nlbf;ii++)
          {
            index = tt[ii];

            bb1 = values1[index];

            F[0] += bb1*dN_dx[ii];
            F[1] += bb1*dN_dy[ii];
          }
          


        for(ii=0;ii<nlbf;ii++)
        {
          bb1 = dN_dx[ii]*dvol0;
          bb2 = dN_dy[ii]*dvol0;
          bb3 = N[ii]*dvol0;

          Flocal[ii]   += (bb3*force - bb1*F[0] - bb2*F[1]) ;

          for(jj=0;jj<nlbf;jj++)
          {
            Klocal(ii, jj)  +=  (bb1*dN_dx[jj] + bb2*dN_dy[jj]) ;
          }
        }
  }//gp1
  }//gp2
  
  //cout << " mmmmmmmmmmmm " << endl;
  //printStiffnessMatrix();
  //printForceVector();
*/
  return 0;
}


int NurbsElem2DHeatTransfer::calcLoadVector()
{
/*
   for(int ii=0;ii<nsize;ii++)
     stiffness_local[ii].zero();

   resi.zero();

   if(tracflag)
   {
    cout << " element number ... = " << elenum << endl;


      double  *gausspoints1 = &(surf0->gausspoints1[0]);
      double  *gausspoints2 = &(surf0->gausspoints2[0]);
      double  *gaussweights1 = &(surf0->gaussweights1[0]);
      double  *gaussweights2 = &(surf0->gaussweights2[0]);
      double  *values1 = &(surf1->Values[0][0]);
   
      int *tt = &(surf0->IEN[elenum][0]);

      double  J, Jmod, b1, b2, params[2], val1, val2, res, dircos[2];
        
      int p = surf0->p, q = surf0->q, ii, jj, gp, TI, TIp1, index, nlbf2;
      double   NN[nlbf], dN_dx[nlbf], dN_dy[nlbf], temp[nlbf], Jac, alpha;

      alpha = elmDat[7];

        // side #1
        if(!CompareDoubles(tracdata[0][0],0.0) )
        {
           nlbf2 = p+1;
           double  N[nlbf2];

           val1 = 0.5*uvalues[2];
           val2 = 0.5*(uvalues[0]+uvalues[1]);

           params[1] = 0.0;
           for(gp=0;gp<nGP1;gp++)   // loop over Gauss points
           {
              params[0] = val1*gausspoints1[gp] + val2;
           
              NurbsShapeFunctions2DAlg2(surf0, startindex[0], startindex[1], gausspoints1[gp], -1.0, N, J, dircos);

              surf0->ShapeFunDerivatives(&(startindex[0]), params, NN, dN_dx, dN_dy, Jac);

              Jmod = J * gaussweights1[gp] ;

              res = 0.0;
              for(ii=0;ii<nlbf;ii++)
              {
                 b1 = values1[tt[ii]];
                 b2 = dN_dy[ii] - alpha*NN[ii];

                 res -= b1 * b2;
                 temp[ii] = b2;
              }
              
              for(ii=0;ii<nsize;ii++)
              {
                b1 = Jmod * temp[ii];
                resi[ii] += b1*res;

                for(jj=0;jj<nsize;jj++)
                  stiffness_local[ii][jj]  +=  b1*temp[jj] ;
              }
           }
           //tracdata[0][0] = 1.0;
           cout << " side1 done " << endl;
        }
        // side #3
        if(!CompareDoubles(tracdata[2][0],0.0) )
        {
           //tracdata[2][0] = 0.0;

           nlbf2 = p+1;
           double   N[nlbf2];

           val1 = 0.5*uvalues[2];
           val2 = 0.5*(uvalues[0]+uvalues[1]);

           params[1] = 1.0;

           for(gp=0;gp<nGP1;gp++)   // loop over Gauss points
           {
              params[0] = val1*gausspoints1[gp] + val2;
           
              NurbsShapeFunctions2DAlg2(surf0, startindex[0], startindex[1], gausspoints1[gp], 1.0, N, J, dircos);

              surf0->ShapeFunDerivatives(&(startindex[0]), params, NN, dN_dx, dN_dy, Jac);

              Jmod = J * gaussweights1[gp] ;
              //Jmod = Jac * gaussweights1[gp] * val1;
              //Jmod = Jac * gaussweights1[gp] * JacMultFact;

              res = -20.0;
              for(ii=0;ii<nlbf;ii++)
              {
                 b1 = values1[tt[ii]];
                 b2 = dN_dy[ii] - alpha*NN[ii];

                 res -= b1 * b2;
                 temp[ii] = b2;
              }
              
              for(ii=0;ii<nsize;ii++)
              {
                b1 = Jmod * temp[ii];
                resi[ii] += b1*res;

                for(jj=0;jj<nsize;jj++)
                  stiffness_local[ii][jj]  +=  b1*temp[jj] ;
              }
           }
           //tracdata[2][0] = 1.0;
           cout << " side3 done " << endl;
        }
        // side #2
        if(!CompareDoubles(tracdata[1][0],0.0) )
        {
           //tracdata[1][0] = 0.0;

           nlbf2 = q+1;
           double   N[nlbf2];

           val1 = 0.5*vvalues[2];
           val2 = 0.5*(vvalues[0]+vvalues[1]);

           params[0] = 1.0;

           for(gp=0;gp<nGP2;gp++)   // loop over Gauss points
           {
              params[1] = val1*gausspoints2[gp] + val2;
           
              NurbsShapeFunctions2DAlg2(surf0, startindex[0], startindex[1], 1.0, gausspoints2[gp], N, J, dircos);

              surf0->ShapeFunDerivatives(&(startindex[0]), params, NN, dN_dx, dN_dy, Jac);

              Jmod = J * gaussweights2[gp] ;
              //Jmod = Jac * gaussweights2[gp] * val1;
              //Jmod = Jac * gaussweights2[gp] * JacMultFact;

              res = 0.0;
              for(ii=0;ii<nlbf;ii++)
              {
                 b1 = values1[tt[ii]];
                 b2 = dN_dx[ii] + alpha*NN[ii];

                 res -= b1 * b2;
                 temp[ii] = b2;
              }
              
              for(ii=0;ii<nsize;ii++)
              {
                b1 = Jmod * temp[ii];
                resi[ii] += b1*res;

                for(jj=0;jj<nsize;jj++)
                  stiffness_local[ii][jj]  +=  b1*temp[jj] ;
              }
           }
           //tracdata[1][0] = 1.0;
           cout << " side2 done " << endl;
        }
        // side #4
        if(!CompareDoubles(tracdata[3][0],0.0) )
        {
           //tracdata[3][0] = 0.0;

           nlbf2 = q+1;
           double   N[nlbf2];

           val1 = 0.5*vvalues[2];
           val2 = 0.5*(vvalues[0]+vvalues[1]);

           params[0] = 0.0;

           for(gp=0;gp<nGP2;gp++)   // loop over Gauss points
           {
              params[1] = val1*gausspoints2[gp] + val2;
           
              NurbsShapeFunctions2DAlg2(surf0, startindex[0], startindex[1], -1.0, gausspoints2[gp], N, J, dircos);

              surf0->ShapeFunDerivatives(&(startindex[0]), params, NN, dN_dx, dN_dy, Jac);

              Jmod = J * gaussweights2[gp] ;
              //Jmod = Jac * gaussweights2[gp] * val1;
              //Jmod = Jac * gaussweights2[gp] * JacMultFact;

              res = 0.0;
              for(ii=0;ii<nlbf;ii++)
              {
                 b1 = values1[tt[ii]];
                 b2 = dN_dx[ii] ;

                 res -= b1 * b2;
                 temp[ii] = b2;
              }
              
              for(ii=0;ii<nsize;ii++)
              {
                b1 = Jmod * temp[ii];
                resi[ii] += b1*res;

                for(jj=0;jj<nsize;jj++)
                  stiffness_local[ii][jj]  +=  b1*temp[jj] ;
              }
           }
           //tracdata[3][0] = 1.0;
           cout << " side4 done " << endl;
        }
//printForceVector();
//printf("\n\n");
//printStiffnessMatrix();
    }
*/
  return 0;
}





int NurbsElem2DHeatTransfer::calcInternalForces()
{
  return 0;
}



void NurbsElem2DHeatTransfer::discreteContourplot(int vartype, int varindex, int index, int nCol, double umin, double umax)
{
  if(index > nivGP)
  {
     cout << '\t' << " Error in NurbsElem2DHeatTransfer::contourplot " << endl;
     return;
  }
  
  if(varindex > 4)
    varindex -= 2;

   vector<double>  outval(nGP);

   switch(vartype)
   {
       case 0:  // plot total strain
       case 1:  // plot elastic strain
       case 2:  // plot plastic strain

                //projectStrain(vartype, varindex, outval);

              break;

       case 3:  // plot stress

                //projectStress(varindex, outval);

              break;

       case 4:  // plot element internal variables

                //projectIntVar(index, outval);

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








int NurbsElem2DHeatTransfer::calcStiffnessMatrix(double dt)
{

  return 0;
}



int NurbsElem2DHeatTransfer::calcOutput(double u1, double v1)
{

  return 0;
}




void NurbsElem2DHeatTransfer::toPostprocess(int vartype, int varindex, int type, SparseMatrixXd&  coeffMat, VectorXd& rhsVec)
{
/*
   double F[4], detF=0.0, F33, Jac, dt, dN_dx[nlbf], dN_dy[nlbf], stre[4], cc[4][4];

   int   err,  isw,  count,  count1, index, ll = 0, ii, jj, gp1, gp2, row, col;

   MatrixXd  Nlocal(nlbf,nlbf);
   VectorXd  NN(nlbf), rhslocal(nlbf);

   Nlocal.setZero();   
   rhslocal.setZero();

   double *gaussweights = &(surf0->gaussweights[0]);

   count = 1;   ll = 0;   err = 0;   isw = 3;
   dt = mpapTime.dt;

   count1 = 0;
   for(gp2=0;gp2<nGP2;gp2++)
   {
   for(gp1=0;gp1<nGP1;gp1++)
   {
        index = 2*count1;

        surf0->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), &NN(0), dN_dx, dN_dy, Jac);

        surf1->deformationGradient(startindex[0], startindex[1], 1, dN_dx, dN_dy, F, detF);

        if(sss == 1)  // plane stress
        {
          if(finite)
            F33 = 1.0/sqrt(detF);
          else
            F33 = 3.0 - F[0] - F[3];
        }
        else if(sss == 2)    // plane strain
          F33 = 1.0;

        matlib2d_(matDat, F, &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);
        count++;
        count1++;
        ll += nivGP;
        
        Nlocal += NN * NN.transpose();
        
        if(varindex < 4)
           rhslocal += NN * stre[varindex];
        else if(varindex == 4)
           rhslocal += NN * sqrt((pow(stre[0]-stre[1],2.0) + pow(stre[1]-stre[2], 2.0) + pow(stre[2]-stre[0], 2.0) + 6.0 * stre[3]*stre[3])/2.0);
        else if(varindex == 5)
           rhslocal += NN * (stre[0]+stre[1]+stre[2])/3.0;

  }//gp1
  }//gp2


    int *tt;
    tt = &(surf0->IEN[elenum][0]);
      
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
*/

  return;
}





void NurbsElem2DHeatTransfer::projectToKnots(bool extrapolateFlag, int vartype, int varindex, int index)
{
/*
   vals2project[0] = intVar2[indx];
   vals2project[1] = intVar2[(nGP1-1)*nivGP+indx];
   vals2project[2] = intVar2[nGP1*(nGP2-1)*nivGP+indx];
   vals2project[3] = intVar2[(nGP1*nGP2-1)*nivGP+indx];
*/

  double  outval[500];

   switch(vartype)
   {
       case 1:  // plot total strain
       case 2:  // plot elastic strain
       case 3:  // plot plastic strain

                projectStrain(vartype, varindex, outval);

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







void NurbsElem2DHeatTransfer::projectStrain(int vartype, int varindex, double* outval)
{
/*
   if(varindex > 5)
   {
       cout << '\t' << "    NurbsElem2DHeatTransfer::projectStrain .... : Error in 'varindex' " << endl;
       return;
   }
   if(vartype > 3)
   {
       cout << '\t' << "    NurbsElem2DHeatTransfer::projectStrain .... : Error in 'vartype' " << endl;
       return;
   }

   double  N[nlbf], dN_dx[nlbf], dN_dy[nlbf], Jac, grad;

   int  count1, ii, gp1, gp2, index;

   double  *values1 = &(surf1->Values[0][0]);
   
   int *tt = &(surf0->IEN[elenum][0]);

   count1 = 0;
   for(gp2=0;gp2<nGP2;gp2++)
   {
   for(gp1=0;gp1<nGP1;gp1++)
   {
        index = count1*2;

        surf0->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, Jac);

        grad = 0.0;

        switch(varindex)
        {
            case 0:    // dT/dx

               for(ii=0;ii<nlbf;ii++)
                 grad += values1[tt[ii]] * dN_dx[ii];

            outval[count1] = grad;

            break;

            case 1:    // dT/dy

               for(ii=0;ii<nlbf;ii++)
                 grad += values1[tt[ii]] * dN_dy[ii];

            outval[count1] = grad;

            break;
        }
        count1++;
  }//gp1
  }//gp2
*/
  return;
}



void  NurbsElem2DHeatTransfer::AssembleElementMatrix(int index, SparseMatrixXd& mtx, int start1, int start2)
{
    int  ii, jj, nn=0, aa, bb, ind, *tt1, *tt2;

    tt1 = &(surf0->LM[elenum][0]);

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







int NurbsElem2DHeatTransfer::calcError(int ind)
{
/*
   int ii, jj, gp1, gp2, count1,  TI, index;
   
   double rad, theta, Jac, fact, dvol0;
   
   double  N[nlbf], dN_dx[nlbf], dN_dy[nlbf], v1, v2, disp[2];
   
   LaplaceNew  analy;
   
   EPOINT  EP;

   double *gaussweights1 = &(surf0->gaussweights1[0]);
   double *gaussweights2 = &(surf0->gaussweights2[0]);

   double  *values1 = &(surf1->Values[0][0]);
   
   int *tt = &(surf0->IEN[elenum][0]);


   elemError = 0.0;

   if(ind == 1) // x - displacement
   {
      count1 = 0;
      for(gp2=0;gp2<nGP2;gp2++)
      {
      for(gp1=0;gp1<nGP1;gp1++)
      {
          index = count1*2;

          surf0->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, Jac);
      
          fact = gaussweights2[gp2] * gaussweights1[gp1] * thick * JacMultFact;

          dvol0 = Jac * fact;
          count1++;

          EP = surf0->SurfacePoint(knotsAtGPs[index], knotsAtGPs[index+1]).CalcEuclid();

          v1 = analy.value(EP.x, EP.y);

          v2 = 0.0;
          for(ii=0;ii<nlbf;ii++)
          {
            v2 += values1[tt[ii]]* N[ii];
          }

          //printf(" \t %12.8f \t %12.8f \n ", v1, v2);

          v1 -= v2;
          
          fact = v1*v1;
          
          elemError += ( fact * dvol0 );
      }
      }
   }
   if(ind == 2) // y - displacement
   {
      count1 = 0;
      for(gp2=0;gp2<nGP2;gp2++)
      {
      for(gp1=0;gp1<nGP1;gp1++)
      {
          index = count1*2;

          surf0->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, Jac);
      
          fact = gaussweights2[gp2] * gaussweights1[gp1] * thick * JacMultFact;

          dvol0 = Jac * fact;
          count1++;

          EP = surf0->SurfacePoint(knotsAtGPs[index], knotsAtGPs[index+1]).CalcEuclid();

          //v1 = analy.dispY(rad, theta);

          //printf(" \t %12.8f \t %12.8f \n ", v1, v2);
          
          v1 -= v2;

          fact = v1*v1;
          
          elemError += ( fact * dvol0 );
      }
      }
   }
   if(ind == 3) // stress
   {
   }
   if(ind == 4) // Energy norm
   {

   }
   if(ind == 5) // LS functional
   {

   }

  //if( (int) matDat[4] == 1)
    //printf(" \t element = %5d ... \t ... elemError  =   %12.6E \n " , elenum, elemError);
*/
  return 0;
}





