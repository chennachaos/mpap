
#include <math.h>
#include "Debug.h"
#include "MpapTime.h"
//#include "Plot.h"
#include "NurbsElem2DStructSolid.h"
#include "NurbsShapeFunctions.h"
#include <assert.h>
#include "ComputerTime.h"
#include "TimeFunction.h"

#include "ExactSolutionsElasticity.h"

using namespace std;

extern List<TimeFunction> timeFunction;
extern ComputerTime       computerTime;
extern MpapTime mpapTime;
//extern Plot plot;


NurbsElem2DStructSolid::NurbsElem2DStructSolid()
{
  if (debug) cout << " constructor NurbsElem2DStructSolid\n\n";
}



NurbsElem2DStructSolid::~NurbsElem2DStructSolid()
{
  if (debug) cout << " destructor NurbsElem2DStructSolid\n\n";
}


int NurbsElem2DStructSolid::calcStiffnessAndResidual()
{
  if(axsy) // for axisymmetric case
    return NurbsElem2DStructSolid::calcStiffnessAndResidual2();
  else  // for plane stress/strain cases
    return NurbsElem2DStructSolid::calcStiffnessAndResidual1();
}



int NurbsElem2DStructSolid::calcInternalForces()
{
  if(axsy) // for axisymmetric case
    return NurbsElem2DStructSolid::calcInternalForcesAxsy();
  else  // for plane stress/strain cases
    return NurbsElem2DStructSolid::calcInternalForces1();
}




int NurbsElem2DStructSolid::calcStiffnessAndResidual1()
{
  double  F[4], detF=0.0, F33, fact, fact1, fact2, dvol, dvol0, Jac, dt, bb1, bb2, bb3;
  double  N[nlbf], dN_dx[nlbf], dN_dy[nlbf], stre[4], cc[4][4], bc[2][4], force[2], bforce[2];
  double  veloCur[2], acceCur[2], acceFact;
  
  EPOINT  EP;

  int   err,  isw,  count,  count1, index, ll = 0, ii, jj, kk;
  int   gp1, gp2, twoI, twoIp1, twoJ, twoJp1;

  //double  BULK,  mu,  E=1.0e5, nu=0.49999;

  //BULK  = E/3.0/(1.0-2.0*nu);
  //mu    = E/2.0/(1.0+nu);

  //matDat[0] = BULK;
  //matDat[1] = mu;
  //matDat[2] = nu;
  
  //ThickCylinder  analy(sss, matDat[1], matDat[2]);
  PlateWithHole  analy(sss, matDat[1], matDat[2]);
  //ElasticityEx2  analy(matDat[0], matDat[1]);
  //ElasticityEx3  analy(matDat[0], matDat[1]);


  Klocal.setZero();
  Flocal.setZero();

  double *gaussweights1 = &(surf0->gaussweights1[0]);
  double *gaussweights2 = &(surf0->gaussweights2[0]);

  double  *valuesCur1 = &(surf1->ValuesCur[0][0]);
  double  *valuesCur2 = &(surf1->ValuesCur[1][0]);

  double  *valuesDotCur1 = &(surf1->ValuesDotCur[0][0]);
  double  *valuesDotCur2 = &(surf1->ValuesDotCur[1][0]);

  double  *valuesDotDotCur1 = &(surf1->ValuesDotDotCur[0][0]);
  double  *valuesDotDotCur2 = &(surf1->ValuesDotDotCur[1][0]);

  double  af = surf1->td(2);
  double  d1 = surf1->td(5);
  double  aa = surf1->td(10);

  double  rho0 = elmDat[7];
  bforce[0] = elmDat[8];
  bforce[1] = elmDat[9];


  int *tt = &(surf0->IEN[elenum][0]);

  count = 1;   ll = 0;   err = 0;   isw = 3;
  dt = mpapTime.dt;

  count1 = 0;
  for(gp2=0;gp2<nGP2;gp2++)
  {
  for(gp1=0;gp1<nGP1;gp1++)
  {
      index = count1*2;

      surf0->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, Jac);

      //for(ii=0; ii<nlbf; ii++)
	//cout << ii << '\t' << dN_dx[ii] << '\t' << dN_dy[ii] << endl;

        fact = gaussweights2[gp2] * gaussweights1[gp1] * thick * JacMultFact;

        dvol0 = Jac * fact;
        dvol = dvol0;

        //EP = surf0->SurfacePoint(knotsAtGPs[index], knotsAtGPs[index+1]).CalcEuclid();

        //force[0] = analy.forceX(EP.x, EP.y);
        //force[1] = analy.forceY(EP.x, EP.y);
        force[0] = 0.0;
        force[1] = 0.0;

        force[0] = rho0*bforce[0]*timeFunction[0].prop;
        force[1] = rho0*bforce[1]*timeFunction[0].prop;

        //cout << " force ... " << force[0] << '\t' << force[1] << endl;

        surf1->deformationGradient(startindex[0], startindex[1], 1, dN_dx, dN_dy, F, detF);

        if(2 < 1)
	{
          F[0] = F[1] = F[2] = F[3] = 0.0;
          for(ii=0;ii<nlbf;ii++)
          {
            kk = tt[ii];

            bb1 = valuesCur1[kk];
            bb2 = valuesCur2[kk];

            F[0] += bb1*dN_dx[ii];
            F[2] += bb1*dN_dy[ii];
            F[1] += bb2*dN_dx[ii];
            F[3] += bb2*dN_dy[ii];
          }
          F[0] += 1.0;
          F[3] += 1.0;
        }
        //detF = F[0]*F[3] - F[1]*F[2];

          veloCur[0] = veloCur[1] = 0.0;
          acceCur[0] = acceCur[1] = 0.0;
          for(ii=0;ii<nlbf;ii++)
          {
            kk = tt[ii];

            veloCur[0] += valuesDotCur1[kk]*N[ii];
            veloCur[1] += valuesDotCur2[kk]*N[ii];

            acceCur[0] += valuesDotDotCur1[kk]*N[ii];
            acceCur[1] += valuesDotDotCur2[kk]*N[ii];
          }

          //cout << " acceCur ... " << acceCur[0] << '\t' << acceCur[1] << endl;

          force[0] -= rho0*acceCur[0];
          force[1] -= rho0*acceCur[1];

        //printf("F... \t%20.18f\t%20.18f\t%20.18f\t%20.18f\t%20.18f \n", F[0], F[1], F[2], F[3], detF);

        if(detF < 0.0)   return 1;

        if(finite)
        {
           surf1->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, Jac);
           dvol = Jac * fact;
        }

      //for(ii=0; ii<nlbf; ii++)
	//cout << ii << '\t' << dN_dx[ii] << '\t' << dN_dy[ii] << endl;

        // ADJUST F33 fOR 2D PROBLEMS BASED ON THE ASSUMPTIONS OF PLANE STRESS/PLANE STRAIN/AXISYMMETRIC

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

        if(err !=0)          return 2;

        if(finite)
          dvol *= F33;

        //printf(" stress... \t%14.12f\t%14.12f\t%14.12f\t%14.12f\n\n", stre[0], stre[1], stre[2], stre[3]);
        //printf(" dvol = %14.12f \n", dvol);
        //
        //for(ii=0;ii<4;ii++)
        //{
          //for(jj=0;jj<4;jj++)            printf("\t %12.8f", cc[ii][jj]);
          //printf("\n");
        //}

        for(ii=0;ii<4;ii++)
        {
          stre[ii] *= dvol;
          for(jj=0;jj<4;jj++)
            cc[ii][jj] *= dvol;
        }

        //==============================================
        // CALCULATE TANGENT STIFFNESS
        //==============================================

        //   part 1. -- material part (not necessarily symmetric!!)

        for(ii=0;ii<nlbf;ii++)
        {
           bb1 = dN_dx[ii];
           bb2 = dN_dy[ii];
           bb3 = N[ii]*dvol;

           bc[0][0] = (bb1 * cc[0][0] + bb2 * cc[3][0]);
           bc[0][1] = (bb1 * cc[0][1] + bb2 * cc[3][1]);
           bc[0][2] = (bb1 * cc[0][3] + bb2 * cc[3][3]);

           bc[1][0] = (bb2 * cc[1][0] + bb1 * cc[3][0]);
           bc[1][1] = (bb2 * cc[1][1] + bb1 * cc[3][1]);
           bc[1][2] = (bb2 * cc[1][3] + bb1 * cc[3][3]);

           twoI   = 2*ii;
           twoIp1 = twoI+1;

           Flocal[twoI]   += ( bb3*force[0] - bb1*stre[0] - bb2*stre[3]) ;
           Flocal[twoIp1] += ( bb3*force[1] - bb1*stre[3] - bb2*stre[1]) ;

           acceFact = d1*rho0*bb3;

           for(jj=0;jj<nlbf;jj++)
           {
              bb1 = dN_dx[jj];
              bb2 = dN_dy[jj];

              twoJ   = 2*jj;
              twoJp1 = twoJ+1;

              Klocal(twoI,  twoJ)    +=  af*(bc[0][0] * bb1 + bc[0][2] * bb2) ;
              Klocal(twoI,  twoJp1)  +=  af*(bc[0][1] * bb2 + bc[0][2] * bb1) ;
              Klocal(twoIp1,twoJ)    +=  af*(bc[1][0] * bb1 + bc[1][2] * bb2) ;
              Klocal(twoIp1,twoJp1)  +=  af*(bc[1][1] * bb2 + bc[1][2] * bb1) ;

              fact = acceFact*N[jj];

              Klocal(twoI,  twoJ)    += fact ;
              Klocal(twoIp1,twoJp1)  += fact ;
           }
        }

        //   part 2. -- geometrical matrix  (if geometry nonlinear)

        if(finite)
        {
          //cout << " aaaaaaaaaaaa " << endl;
          for(ii=0;ii<nlbf;ii++)
           {
              fact1 = dN_dx[ii] * stre[0] + dN_dy[ii] * stre[3] ;
              fact2 = dN_dx[ii] * stre[3] + dN_dy[ii] * stre[1] ;

              twoI   = 2*ii;
              twoIp1 = twoI+1;

              for(jj=0;jj<nlbf;jj++)
              {
                 twoJ   = 2*jj;
                 twoJp1 = twoJ+1;

                 fact = af*(fact1 * dN_dx[jj] + fact2 * dN_dy[jj]);

                 Klocal(twoI,  twoJ)    += fact ;
                 Klocal(twoIp1,twoJp1)  += fact ;
              }
          }
          //cout << " mmmmmmmmmmmm " << endl;
        }
  }//gp1
  }//gp2
  
  //cout << " mmmmmmmmmmmm " << endl;
  //printStiffnessMatrix();
  //printForceVector();

  return 0;
}






int NurbsElem2DStructSolid::calcInternalForces1()
{
/*
   double F[4], detF=0.0, F33, fact, fact1, fact2, dvol, dvol0, Jac, dt, bb1, bb2;

   double  N[nlbf], dN_dx[nlbf], dN_dy[nlbf], stre[4], cc[4][4];

   int   err,  isw,  count,  count1, index, ll = 0, ii, jj, gp1, gp2, twoI;

  Flocal.setZero();

   double *gaussweights = &(surf0->gaussweights[0]);

   count = 1;   ll = 0;   err = 0;   isw = 3;
   dt = mpapTime.dt;

   count1 = 0;
   for(gp2=0;gp2<nGP2;gp2++)
   {
   for(gp1=0;gp1<nGP1;gp1++)
   {
        index = count1*2;

        surf0->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, Jac);

        fact = gaussweights[count1] * thick * JacMultFact;

        dvol0 = Jac * fact;
        dvol = dvol0;

        surf1->deformationGradient(startindex[0], startindex[1], 1, dN_dx, dN_dy, F, detF);

        if(detF < 0.0)   return 1;

        if(finite)
        {
           surf1->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, Jac);
           dvol = Jac * fact;
        }

        // ADJUST F33 fOR 2D PROBLEMS BASED ON THE ASSUMPTIONS OF PLANE STRESS/PLANE STRAIN/AXISYMMETRIC

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

        if(finite)
          dvol *= F33;

        for(ii=0;ii<4;ii++)
          stre[ii] *= dvol;

        for(ii=0;ii<nlbf;ii++)
        {
          bb1 = dN_dx[ii];
          bb2 = dN_dy[ii];

          twoI   = 2*ii;

          Flocal[twoI]   -= (bb1*stre[0] + bb2*stre[3]) ;
          Flocal[twoI+1] -= (bb1*stre[3] + bb2*stre[1]) ;
        }

  }//gp1
  }//gp2
*/

  return 0;
}





int NurbsElem2DStructSolid::calcStiffnessAndResidual2()
{
/*
        //   for axisymmetric problems

   double F[4], detF=0.0, F33=0.0, fact, fact1, fact2, fact3, fact4, stre[4], cc[4][4], bc[2][4];

   double rad0, rad, dvol, dvol0, r1drad, r1drad0, Jac, dt, totvol=0.0, bb1, bb2;

   int   err,  isw,  count,  count1, index, ll = 0, ii, jj, gp1, gp2, twoI, twoIp1, twoJ, twoJp1;

   double N[nlbf], dN_dx[nlbf], dN_dy[nlbf];

   for(ii=0;ii<nsize;ii++)
     stiffness_local[ii].zero();

   resi.zero();

   double *gaussweights = &(surf0->gaussweights[0]);

   count = 1;   ll = 0;   err = 0;   isw = 3;
   dt = mpapTime.dt;

   count1 = 0;
   for(gp2=0;gp2<nGP2;gp2++)
   {
   for(gp1=0;gp1<nGP1;gp1++)
   {
        index = count1*2;

        NurbsShapeFunctions2DAlg11(surf0, startindex[0], startindex[1], knotsAtGPs[index], knotsAtGPs[index+1], N, dN_dx, dN_dy, Jac);

        fact = gaussweights[count1];// JacMultFact;
        dvol0 = Jac * fact;

        rad0 = surf0->SurfacePoint(knotsAtGPs[index], knotsAtGPs[index+1]).CalcEuclid().x;
        rad  = surf1->SurfacePoint(knotsAtGPs[index], knotsAtGPs[index+1]).CalcEuclid().x;

        dvol0 *= (twoPI * rad0);
        dvol = dvol0;
        r1drad0 = 1.0/rad0;
        r1drad  = 1.0/rad;

        //cout << rad0 << '\t' << rad << '\t' << dvol0 << '\t' << dvol << '\t' << Jac << '\t' << thick << endl;

        surf1->deformationGradient(startindex[0], startindex[1], 1, dN_dx, dN_dy, F, detF);

        if(detF < 0.0)   return 1;

        if(finite)
        {
           NurbsShapeFunctions2DAlg11(surf1, startindex[0], startindex[1], knotsAtGPs[index], knotsAtGPs[index+1], N, dN_dx, dN_dy, Jac);
           dvol = Jac * fact * (twoPI * rad);

           for(ii=0;ii<nlbf;ii++)
             N[ii] *= r1drad;
        }
        else
        {
           for(ii=0;ii<nlbf;ii++)
             N[ii] *= r1drad0;
        }

        F33 = rad/rad0;

        matlib2d_(matDat, F, &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);
        count++;
        count1++;
        ll += nivGP;

//        printf(" stresses ");        printf("\t%12.8f\t%12.8f\t%12.8f\t%12.8f\n\n", stre[0], stre[1], stre[2], stre[3]);

          //if(finite)            dvol *= F33;

//          printf(" detF ");        printf("\t%12.8f\t%12.8f\t%12.8f\t%12.8f\n", detF,dvol0,dvol, dvol0*detF);

        for(ii=0;ii<4;ii++)
        {
          stre[ii] *= dvol;
          for(jj=0;jj<4;jj++)
            cc[ii][jj] *= dvol;
        }

        //==============================================
        // CALCULATE TANGENT STIFFNESS
        //==============================================

        //   part 1. -- material part (not necessarily symmetric!!)


        for(ii=0;ii<nlbf;ii++)
        {
            bb1  = dN_dx[ii];
            bb2  = dN_dy[ii];
            fact = N[ii];

            bc[0][0] = (bb1 * cc[0][0] + bb2 * cc[3][0] + fact * cc[2][0]) ;
            bc[0][1] = (bb1 * cc[0][1] + bb2 * cc[3][1] + fact * cc[2][1]);
            bc[0][2] = (bb1 * cc[0][2] + bb2 * cc[3][2] + fact * cc[2][2]);
            bc[0][3] = (bb1 * cc[0][3] + bb2 * cc[3][3] + fact * cc[2][3]);

            bc[1][0] = (bb2 * cc[1][0] + bb1 * cc[3][0]);
            bc[1][1] = (bb2 * cc[1][1] + bb1 * cc[3][1]);
            bc[1][2] = (bb2 * cc[1][2] + bb1 * cc[3][2]);
            bc[1][3] = (bb2 * cc[1][3] + bb1 * cc[3][3]);

            twoI   = 2*ii;
            twoIp1 = twoI+1;

            resi[twoI]   -= (bb1*stre[0] + bb2*stre[3] + fact*stre[2]) ;
            resi[twoIp1] -= (bb1*stre[3] + bb2*stre[1]) ;

            for(jj=ii;jj<nlbf;jj++)
            {
               twoJ   = 2*jj;
               twoJp1 = twoJ+1;

               bb1 = dN_dx[jj];
               bb2 = dN_dy[jj];
               fact = N[jj];

               stiffness_local[twoI][twoJ]     +=  (bc[0][0] * bb1 + bc[0][3] * bb2 + bc[0][2] * fact) ;
               stiffness_local[twoI][twoJp1]   +=  (bc[0][1] * bb2 + bc[0][3] * bb1) ;
               stiffness_local[twoIp1][twoJ]   +=  (bc[1][0] * bb1 + bc[1][3] * bb2 + bc[1][2] * fact) ;
               stiffness_local[twoIp1][twoJp1] +=  (bc[1][1] * bb2 + bc[1][3] * bb1) ;
            }
        }

        //   part 2. -- geometrical matrix  (if geometry nonlinear)

        if(finite)
        {
           for(ii=0;ii<nlbf;ii++)
           {
              fact1 = (dN_dx[ii] * stre[0] + dN_dy[ii] * stre[3]) ;
              fact2 = (dN_dx[ii] * stre[3] + dN_dy[ii] * stre[1]) ;

              fact4 = N[ii]*stre[2];

              twoI   = 2*ii;
              twoIp1 = twoI+1;

              for(jj=ii;jj<nlbf;jj++)
              {
                 twoJ   = 2*jj;
                 twoJp1 = twoJ+1;

                 fact = fact1 * dN_dx[jj] + fact2 * dN_dy[jj];

                 stiffness_local[twoI][twoJ]     += (fact + fact4 * N[jj]);
                 stiffness_local[twoIp1][twoJp1] += fact ;
              }
           }
        }


  }//gp1
  }//gp2


  for(ii=0;ii<nsize;ii++)
  {
     for(jj=ii+1;jj<nsize;jj++)
     {
        stiffness_local[jj][ii] = stiffness_local[ii][jj];
     }
  }

//  printStiffnessMatrix();

cout << '\t' << " Total Volume for element # " << elenum << " is = " << totvol << endl; cout << endl;

//   computerTime.stopAndPrint(fct);
*/
  return 0;
}




int NurbsElem2DStructSolid::calcInternalForcesAxsy()
{
/*
   //   for axisymmetric problems

   double F[4], detF=0.0, F33=0.0, fact, fact1, fact2, fact3, fact4, stre[4], cc[4][4];

   double rad0, rad, dvol, dvol0, r1drad, r1drad0, Jac, dt, totvol=0.0, bb1, bb2;

   int   err,  isw,  count,  count1, index, ll = 0, ii, jj, gp1, gp2, twoI, twoIp1;

   double N[nlbf], dN_dx[nlbf], dN_dy[nlbf];

   resi.zero();

   double *gaussweights = &(surf0->gaussweights[0]);

   count = 1;   ll = 0;   err = 0;   isw = 3;
   dt = mpapTime.dt;

   count1 = 0;
   for(gp2=0;gp2<nGP2;gp2++)
   {
   for(gp1=0;gp1<nGP1;gp1++)
   {
        index = count1*2;

        NurbsShapeFunctions2DAlg11(surf0, startindex[0], startindex[1], knotsAtGPs[index], knotsAtGPs[index+1], N, dN_dx, dN_dy, Jac);

        fact = gaussweights[count1];// JacMultFact;
        dvol0 = Jac * fact;

        rad0 = surf0->SurfacePoint(knotsAtGPs[index], knotsAtGPs[index+1]).CalcEuclid().x;
        rad  = surf1->SurfacePoint(knotsAtGPs[index], knotsAtGPs[index+1]).CalcEuclid().x;

        dvol0 *= (twoPI * rad0);
        dvol = dvol0;
        r1drad0 = 1.0/rad0;
        r1drad  = 1.0/rad;

        //cout << rad0 << '\t' << rad << '\t' << dvol0 << '\t' << dvol << '\t' << Jac << '\t' << thick << endl;

        surf1->deformationGradient(startindex[0], startindex[1], 1, dN_dx, dN_dy, F, detF);

        if(detF < 0.0)   return 1;

        if(finite)
        {
           NurbsShapeFunctions2DAlg11(surf1, startindex[0], startindex[1], knotsAtGPs[index], knotsAtGPs[index+1], N, dN_dx, dN_dy, Jac);
           dvol = Jac * fact * (twoPI * rad);

           for(ii=0;ii<nlbf;ii++)
             N[ii] *= r1drad;
        }
        else
        {
           for(ii=0;ii<nlbf;ii++)
             N[ii] *= r1drad0;
        }

        F33 = rad/rad0;

        matlib2d_(matDat, F, &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);
        count++;
        count1++;
        ll += nivGP;

//        printf(" stresses ");        printf("\t%12.8f\t%12.8f\t%12.8f\t%12.8f\n\n", stre[0], stre[1], stre[2], stre[3]);

          //if(finite)            dvol *= F33;

//          printf(" detF ");        printf("\t%12.8f\t%12.8f\t%12.8f\t%12.8f\n", detF,dvol0,dvol, dvol0*detF);

        for(ii=0;ii<4;ii++)
        {
          stre[ii] *= dvol;
        }

        //==============================================
        // CALCULATE TANGENT STIFFNESS
        //==============================================

        //   part 1. -- material part (not necessarily symmetric!!)


        for(ii=0;ii<nlbf;ii++)
        {
            bb1  = dN_dx[ii];
            bb2  = dN_dy[ii];

            twoI   = 2*ii;
            twoIp1 = twoI+1;

            resi[twoI]   -= (bb1*stre[0] + bb2*stre[3] + N[ii]*stre[2]) ;
            resi[twoIp1] -= (bb1*stre[3] + bb2*stre[1]) ;
        }

  }//gp1
  }//gp2
*/
  return 0;
}





void NurbsElem2DStructSolid::discreteContourplot(int vartype, int varindex, int index, int nCol, double umin, double umax)
{
  if(index > nivGP)
  {
     cout << '\t' << " Error in NurbsElem2DStructSolid::contourplot " << endl;
     return;
  }
  
  if(varindex > 4)
    varindex -= 2;

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
          //plot.triangleContourPlot(x1, x2, x3, u1, u1, u1, umin, umax, nCol);

          // contour plot for 2nd triangle
          //plot.triangleContourPlot(x1, x3, x4, u1, u1, u1, umin, umax, nCol);

          count++;
      }
  }

  return;
}




void NurbsElem2DStructSolid::projectToKnots(bool extrapolateFlag, int vartype, int varindex, int index)
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




void NurbsElem2DStructSolid::projectStress(int varindex, double* outval)
{
/*
   if(varindex > 5)
   {
       cout << '\t' << "    NurbsElem2DStructSolid::projectStress .... : Error in 'varindex' " << endl;
       return;
   }

   double F[4], detF=0.0, F33=0.0, stre[4], cc[4][4], N[nlbf], dN_dx[nlbf], dN_dy[nlbf];

   double rad0=0.0, rad=0.0, Jac, dt;

   int  err, isw, count, count1, ll, ii, jj, index, gp1, gp2;

   int nivEL = nGP * nivGP;
   for(ii=0;ii<nivEL;ii++)
     intVar2[ii] = intVar1[ii];

   count = 1;   ll = 0;   err = 0;   isw = 3;
   dt = mpapTime.dt;

   count1 = 0;
   for(gp2=0;gp2<nGP2;gp2++)
   {
   for(gp1=0;gp1<nGP1;gp1++)
   {
        index = count1*2;

        surf0->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, Jac);

        surf1->deformationGradient(startindex[0], startindex[1], 1, dN_dx, dN_dy, F, detF);

        //   for axisymmetric problems compute radius
        if(axsy)
        {
           rad0 = surf0->SurfacePoint(knotsAtGPs[index], knotsAtGPs[index+1]).CalcEuclid().x;
           rad  = surf1->SurfacePoint(knotsAtGPs[index], knotsAtGPs[index+1]).CalcEuclid().x;
        }

        if(sss == 1)  // plane stress
        {
          if(finite)
            F33 = 1.0/sqrt(detF);
          else
            F33 = 3.0 - F[0] - F[3];
        }
        else if(sss == 2)    // plane strain
          F33 = 1.0;
        else //(sss == 3)    // axisymmetric
          F33 = rad/rad0;

        matlib2d_(matDat, F, &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, (Element*) this);
        
        //printf("\t%12.8f\t%12.8f\t%12.8f\n",stre[0], stre[1], stre[2]);

        if(varindex < 4)
           outval[count1] = stre[varindex];
        else if(varindex == 4)
           outval[count1] = sqrt((pow(stre[0]-stre[1],2.0) + pow(stre[1]-stre[2], 2.0) + pow(stre[2]-stre[0], 2.0) + 6* stre[3]*stre[3])/2);
        else if(varindex == 5)
           outval[count1] = (stre[0]+stre[1]+stre[2])/3.0;

        count++;
        count1++;
        ll += nivGP;

  }//gp1
  }//gp2
*/
  return;
}




void NurbsElem2DStructSolid::projectStrain(int vartype, int varindex, double* outval)
{
/*
   if(varindex > 5)
   {
       cout << '\t' << "    NurbsElem2DStructSolid::projectStress .... : Error in 'varindex' " << endl;
       return;
   }
   if(vartype > 3)
   {
       cout << '\t' << "    NurbsElem2DStructSolid::projectStress .... : Error in 'vartype' " << endl;
       return;
   }

   double F[4], detF=0.0, F33=0.0, stre[4], cc[4][4], N[nlbf], dN_dx[nlbf], dN_dy[nlbf];

   double rad0=0.0, rad=0.0, Jac, dt;

   int  err, isw, count, count1, ll, ii, jj, index, gp1, gp2;

   int nivEL = nGP * nivGP;
   for(ii=0;ii<nivEL;ii++)
     intVar2[ii] = intVar1[ii];

   bool intVarFlag = false;
   if( nivGP > 0)
      intVarFlag = true;

   count = 1;   ll = 0;   err = 0;   isw = 3;
   dt = mpapTime.dt;

   count1 = 0;
   for(gp2=0;gp2<nGP2;gp2++)
   {
   for(gp1=0;gp1<nGP1;gp1++)
   {
        index = count1*2;

        surf0->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, Jac);

        surf1->deformationGradient(startindex[0], startindex[1], 1, dN_dx, dN_dy, F, detF);

        if(detF < 0.0)
        {
           cerr << "        NurbsElem2DStructSolid::projectStress.......Negative DEFORMATION GRADIENT   " << endl;
           return ;
        }

        //   for axisymmetric problems compute radius
        if(axsy)
        {
           rad0 = surf0->SurfacePoint(knotsAtGPs[index], knotsAtGPs[index+1]).CalcEuclid().x;
           rad  = surf1->SurfacePoint(knotsAtGPs[index], knotsAtGPs[index+1]).CalcEuclid().x;
        }

        if(sss == 1)  // plane stress
        {
          if(finite)
            F33 = 1.0/sqrt(detF);
          else
            F33 = 3.0 - F[0] - F[3];
        }
        else if(sss == 2)    // plane strain
          F33 = 1.0;
        else //(sss == 3)    // axisymmetric
          F33 = rad/rad0;

        matlib2d_(matDat, F, &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, (Element*) this);

        switch(vartype)
        {
             case 1: // total strain components

                     switch(varindex)
                     {
                          case 0:    // epsxx

                              outval[count1] = F[0] - 1.0;           break;

                          case 1:    // epsyy

                              outval[count1] = F[3] - 1.0;           break;

                          case 2:    // epsxy

                              outval[count1] = F[1] + F[2];           break;
                     }

               break;

             case 2: // plastic strain components

                 if(intVarFlag)
                 {
                     switch(varindex)
                     {
                          case 0:    // epspxx

                              outval[count1] = intVar2[0];         break;

                          case 1:    // epspyy

                              outval[count1] = intVar2[1];         break;

                          case 2:    // epspxy

                              outval[count1] = intVar2[3];         break;

                          case 3:    // epspeqv

                              outval[count1] = intVar2[6];         break;
                     }
                 }

               break;

             case 3: // elastic strain components

               if(intVarFlag)
               {
                     switch(varindex)
                     {
                          case 0:    // epspxx

                              outval[count1] = F[0] - 1.0 - intVar2[0];         break;

                          case 1:    // epspyy

                              outval[count1] = F[3] - 1.0 - intVar2[1];         break;

                          case 2:    // epspxy

                              outval[count1] = F[1] + F[2] - intVar2[3];         break;
                     }

               }
               else
               {
                     switch(varindex)
                     {
                          case 0:    // epsxx

                              outval[count1] = F[0] - 1.0;           break;

                          case 1:    // epsyy

                              outval[count1] = F[3] - 1.0;           break;

                          case 2:    // epsxy

                              outval[count1] = F[1] + F[2];           break;
                     }
               }

               break;
        }

        count++;
        count1++;
        ll += nivGP;

  }//gp1
  }//gp2
*/
  return;
}






void NurbsElem2DStructSolid::projectIntVar(int index, double* outval)
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



int NurbsElem2DStructSolid::calcStiffnessMatrix(double dt)
{

  return 0;
}









int NurbsElem2DStructSolid::calcMassMatrix(int lumpInd, double dt)
{
   /*
   *  lumpInd = 1 --> consistent Mass matrix

   *          = 2 --> Row-sum Mass lumping
   *          = 3 --> proportional Mass lumping
   */

/*
   double fact1, fact2, fact, dvol0=0.0, rad0=0.0, J=0.0;

   VectorArray<double> NN;   NN.setDim(nlbf);


    mass_local.setDim(nsize);
    for(int ii=0;ii<nsize;ii++)
    {

       mass_local[ii].setDim(nsize);
       mass_local[ii].zero();

    }


       // loop over Gauss points

       for(int gp1=0;gp1<gausspoints1.n;gp1++)
       {
          for(int gp2=0;gp2<gausspoints2.n;gp2++)
          {

              NN.zero();

              NurbsShapeFunctions2DAlg3(surf0, startindex[0], startindex[1], gausspoints1[gp1], gausspoints2[gp2], NN, J);

              dvol0 = J * gaussweights[count1] * thick;

              //   for axisymmetric problems compute radius
              if(axsy)

              {
                 rad0 = 0.0 ;
                 int loc_num=0;
                 for(int jj=0; jj<=surf0->q; jj++)

                 {
                    for(int ii=0; ii<=surf0->p; ii++)
                    {
                       rad0 += (NN[loc_num] * surf0->Pw[startindex[0]+ii][startindex[1]+jj].CalcEuclid().x);

                       loc_num++;
                    }
                 }
                 dvol0 *= (twoPI * rad0);

              }

              for(int ii=0;ii<nlbf;ii++)
              {

                  fact = NN[ii] * dvol0;
                  for(int jj=0;jj<nlbf;jj++)
                  {
                     fact1 = fact * NN[jj];

                     mass_local[2*ii][2*jj]     += fact1;
                     mass_local[2*ii+1][2*jj+1] += fact1;
                  }
              }

          }
       }
*/
  return 0;
}





int NurbsElem2DStructSolid::calcOutput(double u1, double v1)
{
/*
   u1 = v1 = gausspoints[0];

   double F[4], detF=0.0, F33=0.0, fact, fact1, fact2, fact3, fact4,
          stre[4], cc[4][4], bc[2][4];

   double rad0=0.0, rad=0.0, dvol=0.0, dvol0=0.0, r1dF33=0.0, r1drad0=0.0, dt, J;

   int   err   = 0,
         isw   = 3,
         count = 1,
         ll    = 0,
         nivEl, ii, jj;

   nivEl = nivGP * nGP;

    for(ii=0;ii<nivEl;ii++)
       intVar2[ii] = intVar1[ii];

   double  N[nlbf], dN_dx[nlbf], dN_dy[nlbf];

        NurbsShapeFunctions2DAlg11(surf0, startindex[0], startindex[1], knotsAtGPs[index], knotsAtGPs[index+1], N, dN_dx, dN_dy, J);

        //   for axisymmetric problems compute radius
        if(axsy)
        {
           rad0 = surf0->SurfacePoint(u1, v1).CalcEuclid().x;
           rad  = surf1->SurfacePoint(u1, v1).CalcEuclid().x;

           r1drad0 = 1.0 / rad0;
        }

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
        else //(sss == 3)    // axisymmetric
          F33 = rad/rad0;

        dt = mpapTime.dt;
        matlib2d_(matDat, F, &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, (Element*) this);

        if(err !=0)
          return 1;
      //
         cout << "      stress vector ...:  " ;
         for(int ii=0;ii<4;ii++)
         {
           cout << fixed << '\t' << stre[ii];
         }
         cout << endl;
      //
*/
  return 0;
}



void NurbsElem2DStructSolid::toPostprocess(int vartype, int varindex, int type, SparseMatrixXd&  coeffMat, VectorXd& rhsVec)
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



/*
void  NurbsElem2DStructSolid::AssembleElementMatrix(int index, Mat mtx, int start1, int start2)
{
    PetscErrorCode ierr;
    int  ii, jj, nn=0, aa, bb, size2, ind, *tt1, *tt2;

    tt1 = &(surf0->LM[elenum][0]);
    
    //cout << elenum << '\t' << elenum2 << endl;
    //cout << nsize << '\t' << size2 << endl;
    //cout << surf0->LM[elenum] << endl;

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
      }
    }

  return;
}
*/


void  NurbsElem2DStructSolid::AssembleElementMatrix(int index, SparseMatrixXd& mtx, int start1, int start2)
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





int NurbsElem2DStructSolid::calcError(int ind)
{
/*
   int ii, jj, gp1, gp2, count1,  TI, index;
   
   double rad, theta, Jac, fact, dvol0, er;
   
   double  N[nlbf], dN_dx[nlbf], dN_dy[nlbf], v1, v2, disp1[2], disp2[2];

   PlateWithHole  analy(sss, matDat[1], matDat[2]);
   //ThickCylinder  analy(sss, matDat[1], matDat[2]);
   //ElasticityEx2  analy(matDat[0], matDat[1]);
   //ElasticityEx1  analy(matDat[1]);
   //ElasticityEx3  analy(matDat[0], matDat[1]);

   EPOINT  EP, EP2;

   double *gaussweights1 = &(surf0->gaussweights1[0]);
   double *gaussweights2 = &(surf0->gaussweights2[0]);

   double  *values1 = &(surf1->Values[0][0]);
   double  *values2 = &(surf1->Values[1][0]);
   
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
          
          EP2 = surf1->SurfacePoint(knotsAtGPs[index], knotsAtGPs[index+1]).CalcEuclid();

          rad = sqrt(EP.x*EP.x + EP.y*EP.y) ;

          theta = atan2(EP.y, EP.x);

          disp1[0] = analy.dispX(rad, theta);
          disp1[1] = analy.dispY(rad, theta);

          //disp1[0] = analy.dispX(EP.x, EP.y);
          //disp1[1] = analy.dispY(EP.x, EP.y);

          disp2[0] = disp2[1] = 0.0;
          //for(ii=0;ii<nlbf;ii++)
          //{
            //disp2[0] += values1[tt[ii]]* N[ii];
            //disp2[1] += values2[tt[ii]]* N[ii];
          //}
          disp2[0] = EP2.x - EP.x;
          disp2[1] = EP2.y - EP.y;

          //printf(" \t %12.8f \t %12.8f \n ", EP.x, EP.y);
          //printf(" \t %12.8f \t %12.8f \t %12.8f \n ", v1, v2, (-disp[0]*sin(theta) + disp[1]*cos(theta)));

          disp1[0] -= disp2[0] ;
          disp1[1] -= disp2[1] ;

          fact = disp1[0]*disp1[0] + disp1[1]*disp1[1];

          elemError += ( fact * dvol0 );
      }
      }
   }
   if(ind == 2) // y - displacement
   {
   }
   if(ind == 3) // stress
   {
      int ll, isw, err, count;

      double F[4], detF, F33, stre[4], stre2[4], cc[4][4], bb1, bb2, dt;

      count = 1;   ll = 0;   err = 0;   isw = 3;
      dt = mpapTime.dt;

      count1 = 0;
      for(gp2=0;gp2<nGP2;gp2++)
      {
      for(gp1=0;gp1<nGP1;gp1++)
      {
          index = count1*2;

          surf0->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, Jac);
      
          fact = gaussweights2[gp2] * gaussweights1[gp1] * thick * JacMultFact;

          dvol0 = Jac * fact;

          EP = surf0->SurfacePoint(knotsAtGPs[index], knotsAtGPs[index+1]).CalcEuclid();

          rad = sqrt(EP.x*EP.x + EP.y*EP.y) ;

          theta = atan2(EP.y, EP.x);
          //theta *= -1.0;
          //cout << " theta = " << (theta*180.0/PI) << endl;

          stre[0] = analy.stressXX(rad, theta);
          stre[1] = analy.stressYY(rad, theta);

          stre[2] = 0.0;
          if(sss == 2)
            stre[2] = matDat[2]*(stre[0]+stre[1]);

          stre[3] = analy.stressXY(rad, theta);

          //stre[0] = analy.stressXX(EP.x, EP.y);
          //stre[1] = analy.stressYY(EP.x, EP.y);
          //stre[2] = 0.0;
          //stre[3] = analy.stressXY(EP.x, EP.y);

          surf1->deformationGradient(startindex[0], startindex[1], 1, dN_dx, dN_dy, F, detF);
          
          if(2 < 1)
	  {
          F[0] = F[1] = F[2] = F[3] = 0.0;
          for(ii=0;ii<nlbf;ii++)
          {
            jj = tt[ii];

            bb1 = values1[jj];
            bb2 = values2[jj];

            F[0] += bb1*dN_dx[ii];
            F[2] += bb1*dN_dy[ii];
            F[1] += bb2*dN_dx[ii];
            F[3] += bb2*dN_dy[ii];
          }
          F[0] += 1.0;
          F[3] += 1.0;
          }

        //printf("F... \t%20.18f\t%20.18f\t%20.18f\t%20.18f\t%20.18f \n", F[0], F[1], F[2], F[3], detF);


        // ADJUST F33 fOR 2D PROBLEMS BASED ON THE ASSUMPTIONS OF PLANE STRESS/PLANE STRAIN/AXISYMMETRIC

        if(sss == 1)  // plane stress
        {
          if(finite)
            F33 = 1.0/sqrt(detF);
          else
            F33 = 3.0 - F[0] - F[3];
        }
        else if(sss == 2)    // plane strain
          F33 = 1.0;

        matlib2d_(matDat, F, &F33, stre2, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);
        count++;
        count1++;
        ll += nivGP;

          for(ii=0;ii<4; ii++)
            stre[ii] -= stre2[ii];

        //printf("stresses... \t%20.18f\t%20.18f\t%20.18f\t%20.18f \n", stre[0], stre[1], stre[2], stre[3]);

          fact = stre[0]*stre[0] + stre[1]*stre[1] + stre[3]*stre[3];// + stre[3]*stre[3] ;
          //fact = stre[0]*stre[0] ;

          elemError += ( fact * dvol0 );
      }
      }
   }
   if(ind == 4) // Energy norm
   {
      int ll, isw, err, count;

      double F1[4], F2[4], eps[4], detF, F33, stre[4], cc[4][4], bb1, bb2, dt;

      count = 1;   ll = 0;   err = 0;   isw = 3;
      dt = mpapTime.dt;

      count1 = 0;
      for(gp2=0;gp2<nGP2;gp2++)
      {
      for(gp1=0;gp1<nGP1;gp1++)
      {
          index = count1*2;

          surf0->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, Jac);
      
          fact = gaussweights2[gp2] * gaussweights1[gp1] * thick * JacMultFact;

          dvol0 = Jac * fact;

          EP = surf0->SurfacePoint(knotsAtGPs[index], knotsAtGPs[index+1]).CalcEuclid();

          rad = sqrt(EP.x*EP.x + EP.y*EP.y) ;

          theta = atan2(EP.y, EP.x);

          F1[0] = analy.strainXX(rad, theta) ;
          F1[1] = analy.strainXY(rad, theta) ;
          F1[2] = F1[1];
          F1[3] = analy.strainYY(rad, theta) ;

          //stre[0] = analy.stressXX(EP.x, EP.y);
          //stre[1] = analy.stressYY(EP.x, EP.y);
          //stre[2] = 0.0;
          //stre[3] = analy.stressXY(EP.x, EP.y);

          //surf1->deformationGradient(startindex[0], startindex[1], 1, dN_dx, dN_dy, F, detF);
          
          //
          F2[0] = F2[1] = F2[2] = F2[3] = 0.0;
          for(ii=0;ii<nlbf;ii++)
          {
            jj = tt[ii];

            bb1 = values1[jj];
            bb2 = values2[jj];

            F2[0] += bb1*dN_dx[ii];
            F2[2] += bb1*dN_dy[ii];
            F2[1] += bb2*dN_dx[ii];
            F2[3] += bb2*dN_dy[ii];
          }
          //F2[0] += 1.0;
          //F2[3] += 1.0;

        for(ii=0;ii<4; ii++)
          F1[ii] -= F2[ii];

        F1[0] += 1.0;
        F1[3] += 1.0;


        //printf("F... \t%20.18f\t%20.18f\t%20.18f\t%20.18f\t%20.18f \n", F[0], F[1], F[2], F[3], detF);


        // ADJUST F33 fOR 2D PROBLEMS BASED ON THE ASSUMPTIONS OF PLANE STRESS/PLANE STRAIN/AXISYMMETRIC

        if(sss == 1)  // plane stress
        {
          if(finite)
            F33 = 1.0/sqrt(detF);
          else
            F33 = 3.0 - F1[0] - F1[3];
        }
        else if(sss == 2)    // plane strain
          F33 = 1.0;


        matlib2d_(matDat, F1, &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);
        count++;
        count1++;
        ll += nivGP;

        eps[0] = F1[0]-1.0;
        eps[1] = F1[3]-1.0;
        eps[2] = 0.0;
        eps[3] = F1[1];

        fact = 0.5*(stre[0]*eps[0] + stre[1]*eps[1] + stre[3]*eps[3] );// + stre[3]*stre[3] ;

        elemError += ( fact * dvol0 );
      }
      }
   }
   if(ind == 5) // LS functional
   {

   }

  //if( (int) matDat[4] == 1)
    //printf(" \t element = %5d ... \t ... elemError  =   %12.6E \n " , elenum, elemError);
*/
  return 0;
}



void NurbsElem2DStructSolid::computeBounds(double* val)
{
  EPOINT  EP1, EP2, EP3, EP4;

  EP1 = surf0->SurfacePoint(uvalues[0], vvalues[0]).CalcEuclid();
  EP2 = surf0->SurfacePoint(uvalues[1], vvalues[0]).CalcEuclid();
  EP3 = surf0->SurfacePoint(uvalues[0], vvalues[1]).CalcEuclid();
  EP4 = surf0->SurfacePoint(uvalues[1], vvalues[1]).CalcEuclid();
  
  EP1 = EP1 - EP4;
  EP2 = EP2 - EP3;
  
  val[0] = max(EP1.Norm(), EP2.Norm());

  return;
}








