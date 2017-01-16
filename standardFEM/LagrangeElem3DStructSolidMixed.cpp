

#include "LagrangeElem3DStructSolidMixed.h"

#include "Debug.h"
#include "MpapTime.h"
#include "ComputerTime.h"
#include "GeomDataLagrange.h"
#include "SolutionData.h"
#include "FunctionsMaterial.h"
#include "NurbsShapeFunctions.h"


using namespace std;

extern ComputerTime       computerTime;
extern MpapTime mpapTime;


LagrangeElem3DStructSolidMixed::LagrangeElem3DStructSolidMixed()
{
  if (debug) cout << " constructor LagrangeElem3DStructSolidMixed\n\n";
}

LagrangeElem3DStructSolidMixed::~LagrangeElem3DStructSolidMixed()
{
  if (debug) cout << " destructor LagrangeElem3DStructSolidMixed\n\n";
}


void LagrangeElem3DStructSolidMixed::prepareElemData()
{
  LagrangeElement::prepareElemData();

   int ii, jj;
   
   npElem = (degree+1)*(degree+1);
   nlbf   = npElem;
   nsize  = npElem * ndof;

   // set the element property variables

    elmDat = &(SolnData->ElemProp[elmType].data[0]);
    matDat = &(SolnData->MatlProp[matType].data[0]);

   finiteInt  = (int) elmDat[2] ;
   sss        = (int) elmDat[3] ;
   thick      = elmDat[4] ;
   //bforce[0]  = elmDat[5] ;
   //bforce[1]  = elmDat[6] ;
   //rho0       = elmDat[7] ;
   matId      = SolnData->MatlProp[matType].id + 1;
   finite     = (finiteInt >= 1) ;
   axsy       = (sss == 3);
   followerLoadFlag = (elmDat[8] == 1);

   if(sss != 1) thick = 1.0; // for plane strain and axisymmetric problems


  //Klocal.resize(nsize, nsize);
  //Flocal.resize(nsize);
  
  return;
}




void LagrangeElem3DStructSolidMixed::prepareElemData2()
{
  return;
}



int LagrangeElem3DStructSolidMixed::calcLoadVector()
{
  return 0;
}


void LagrangeElem3DStructSolidMixed::resetMatrixAndVector()
{
  return;
}


int LagrangeElem3DStructSolidMixed::calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal)
{
  if(finite)
    calcStiffnessAndResidualFS(Klocal, Flocal);
  else
    calcStiffnessAndResidualSS(Klocal, Flocal);

  return 0;
}



int LagrangeElem3DStructSolidMixed::calcStiffnessAndResidualSS(MatrixXd& Klocal, VectorXd& Flocal)
{
  // subroutine for small strain formulation

  //   char fct[] = "LagrangeElem3DStructSolidMixed::calcStiffnessAndResidual2";
  //   computerTime.go(fct);

   int   err,  isw,  count,  count1, index, ll = 0, ii, jj, gp1, gp2, TI, TIp1, TJ, TJp1, sizep, mm;

   double F[4], detF=0.0, F33, fact, fact1, fact2, dvol, dvol0, Jac, dt, totvol=0.0, bb1, bb2, JacMult;

   vector<double>  N(nlbf), dN_dx(nlbf), dN_dy(nlbf), Nbar(sizep);

   double  cc[4][4], bc[2][4], Idev[4][4], D11[4][4], cctmp[4][4], stre[4];
   double  pres, pbar, dummy, strdev[4], BULK, volstr, r2d3=2.0/3.0;
   
   Idev2D(Idev);

   double  *gausspoints1  = &(GeomData->gausspoints1[0]);
   double  *gaussweights1 = &(GeomData->gaussweights1[0]);
   double  *gausspoints2  = &(GeomData->gausspoints2[0]);
   double  *gaussweights2 = &(GeomData->gaussweights2[0]);

   count = 1;   ll = 0;   err = 0;   isw = 3;
   dt = mpapTime.dt;

    if(Klocal.rows() != nsize)
    {
      Klocal.resize(nsize, nsize);
      Flocal.resize(nsize);
    }
    Klocal.setZero();
    Flocal.setZero();


   for(gp2=0;gp2<nGP;gp2++)
   {
     JacMult = gaussweights2[gp2] * thick;

     for(gp1=0;gp1<nGP;gp1++)
     {
        //GeomData->computeBasisFunctions2D(0, gausspoints1[gp1], gausspoints2[gp2], &N[0], &dN_dx[0], &dN_dy[0], Jac);

        fact = gaussweights1[gp1] * JacMult;

        dvol0 = Jac * fact;
        dvol = dvol0;

        //GeomData->computeDeformationGradient(1, dN_dx, dN_dy, F, detF);

        if(detF < 0.0)   return 1;
        
        volstr = F[0] + F[3] - 2.0;

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
        ll += nivGP;

        if(err !=0) 
          return 1;

          //pres = surf2->computeValueAndShanpeFns(1, knotsAtGPs[index], knotsAtGPs[index+1], Nbar);

          dummy = pres - (stre[0]+stre[1]+stre[2])/3.0 ;

          stre[0] += dummy;
          stre[1] += dummy;
          stre[2] += dummy;

//        printf(" stresses ");        printf("\t%12.8f\t%12.8f\t%12.8f\t%12.8f\n\n", stre[0], stre[1], stre[2], pres);

          for(ii=0;ii<4;ii++)
          {
             for(jj=0;jj<4;jj++)
             {
                cctmp[ii][jj] = 0.0;
                for(mm=0;mm<4;mm++)
                   cctmp[ii][jj] += Idev[ii][mm]*cc[mm][jj];
             }
          }

          for(ii=0;ii<4;ii++)
          {
             for(jj=0;jj<4;jj++)
             {
               D11[ii][jj] = 0.0;
               for(mm=0;mm<4;mm++)
                  D11[ii][jj] += cctmp[ii][mm]*Idev[mm][jj];
             }
          }

          for(ii=0;ii<4;ii++)
          {
             stre[ii] *= dvol0;
             for(jj=0;jj<4;jj++)
               D11[ii][jj] *= dvol0;
          }

          //==============================================
          // CALCULATE TANGENT STIFFNESS and RESIDUAL
          //==============================================

          for(ii=0;ii<nlbf;ii++)
          {
             bb1 = dN_dx[ii];
             bb2 = dN_dy[ii];

             bc[0][0] = (bb1 * D11[0][0] + bb2 * D11[3][0]);
             bc[0][1] = (bb1 * D11[0][1] + bb2 * D11[3][1]);
             bc[0][2] = (bb1 * D11[0][3] + bb2 * D11[3][3]);

             bc[1][0] = (bb2 * D11[1][0] + bb1 * D11[3][0]);
             bc[1][1] = (bb2 * D11[1][1] + bb1 * D11[3][1]);
             bc[1][2] = (bb2 * D11[1][3] + bb1 * D11[3][3]);

             TI   = 2*ii;
             TIp1 = TI+1;

             Flocal(TI)   -= (bb1*stre[0] + bb2*stre[3]) ;
             Flocal(TIp1) -= (bb1*stre[3] + bb2*stre[1]) ;

             for(jj=0;jj<nlbf;jj++)
             {
                bb1 = dN_dx[jj];
                bb2 = dN_dy[jj];

                TJ   = 2*jj;
                TJp1 = TJ+1;

                Klocal(TI,TJ)     +=  (bc[0][0] * bb1 + bc[0][2] * bb2) ;
                Klocal(TI,TJp1)   +=  (bc[0][1] * bb2 + bc[0][2] * bb1) ;
                Klocal(TIp1,TJ)   +=  (bc[1][0] * bb1 + bc[1][2] * bb2) ;
                Klocal(TIp1,TJp1) +=  (bc[1][1] * bb2 + bc[1][2] * bb1) ;
             }

             // compute Kup and Kpp matrices

             bb1 *= dvol0;
             bb2 *= dvol0;

             for(jj=0;jj<sizep;jj++)
             {
               Kup(TI, jj)   += ( bb1 * Nbar[jj] );
               Kup(TIp1, jj) += ( bb2 * Nbar[jj] );
             }
          }

          fact1 = (volstr - pres/BULK)*dvol0;
          fact2 = dvol0/BULK;

          for(ii=0;ii<sizep;ii++)
          {
            //Flocal2(ii) -= fact1 * Nbar[ii];

            fact = Nbar[ii] * fact2;

            for(jj=0;jj<sizep;jj++)
               Kpp(ii,jj) -= ( fact * Nbar[jj] );
          }
  }//gp1
  }//gp2

  return 0;
}





int LagrangeElem3DStructSolidMixed::calcStiffnessAndResidualFS(MatrixXd& Klocal, VectorXd& Flocal)
{
  // subroutine for finite strain formulation

  //   char fct[] = "LagrangeElem3DStructSolidMixed::calcStiffnessAndResidual2";
  //   computerTime.go(fct);

   int   err,  isw,  count,  count1, index, ll = 0, ii, jj, gp1, gp2, TI, TIp1, TJ, TJp1, sizep, mm;

   double F[4], detF=0.0, F33, fact, fact1, fact2, fact3, fact4, dvol, dvol0, Jac, dt, totvol=0.0, bb1, bb2, JacMult;

   vector<double>  N(nlbf), dN_dx(nlbf), dN_dy(nlbf), Nbar(sizep);

   double  cc[4][4], bc[2][4], Idev[4][4], D11[4][4], cctmp[4][4], stre[4];
   double  pres, pbar, dummy, strdev[4], BULK, volstr, r2d3=2.0/3.0;

   Idev2D(Idev);

   double  *gausspoints1  = &(GeomData->gausspoints1[0]);
   double  *gaussweights1 = &(GeomData->gaussweights1[0]);
   double  *gausspoints2  = &(GeomData->gausspoints2[0]);
   double  *gaussweights2 = &(GeomData->gaussweights2[0]);

   count = 1;   ll = 0;   err = 0;   isw = 3;
   dt = mpapTime.dt;

    if(Klocal.rows() != nsize)
    {
      Klocal.resize(nsize, nsize);
      Flocal.resize(nsize);
    }
    Klocal.setZero();
    Flocal.setZero();


   for(gp2=0;gp2<nGP;gp2++)
   {
     JacMult = gaussweights2[gp2] * thick;

     for(gp1=0;gp1<nGP;gp1++)
     {
        //GeomData->computeBasisFunctions2D(0, gausspoints1[gp1], gausspoints2[gp2], N, dN_dx, dN_dy, Jac);

        fact = gaussweights1[gp1] * JacMult;

        dvol0 = Jac * fact;
        dvol = dvol0;

        //GeomData->computeDeformationGradient(1, dN_dx, dN_dy, F, detF);

        if(detF < 0.0)
          return 1;

        //GeomData->computeBasisFunctions2D(1, gausspoints1[gp1], gausspoints2[gp2], N, dN_dx, dN_dy, Jac);
        dvol = Jac * fact;

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
        ll += nivGP;

        if(err !=0)    return 1;

        dvol *= F33; // finite strain

          //pres = surf2->computeValueAndShanpeFns(1, knotsAtGPs[index], knotsAtGPs[index+1], Nbar);

          pbar = (stre[0]+stre[1]+stre[2])/3.0;

          strdev[0] = stre[0] - pbar;
          strdev[1] = stre[1] - pbar;
          strdev[2] = stre[2] - pbar;
          strdev[3] = stre[3];

          stre[0] = strdev[0] + pres;
          stre[1] = strdev[1] + pres;
          stre[2] = strdev[2] + pres;

//        printf(" stresses ");        printf("\t%12.8f\t%12.8f\t%12.8f\t%12.8f\n\n", stre[0], stre[1], stre[2], pres);

          for(ii=0;ii<4;ii++)
          {
             for(jj=0;jj<4;jj++)
             {
                cctmp[ii][jj] = 0.0;
                for(mm=0;mm<4;mm++)
                  cctmp[ii][jj] += Idev[ii][mm] * D11[mm][jj];
             }
          }

          for(ii=0;ii<4;ii++)
          {
             for(jj=0;jj<4;jj++)
             {
                D11[ii][jj] = 0.0;
                for(mm=0;mm<4;mm++)
                  D11[ii][jj] += cctmp[ii][mm] * Idev[mm][jj];
             }
          }

          fact = 2.0 * (pbar - pres);
          fact1 = (r2d3*pbar - pres);

          for(ii=0;ii<3;ii++)
          {
             for(jj=0;jj<4;jj++)
             {
                D11[ii][jj] -= r2d3*strdev[jj] ;
                D11[jj][ii] -= r2d3*strdev[jj] ;
             }
             for(jj=0;jj<3;jj++)
                D11[ii][jj] -= fact1;

             D11[ii][ii] += fact;
          }

          D11[3][3] += 0.5*fact;

          for(ii=0;ii<4;ii++)
          {
             stre[ii] *= dvol;
             for(jj=0;jj<4;jj++)
               D11[ii][jj] *= dvol;
          }

          //==============================================
          // CALCULATE TANGENT STIFFNESS and RESIDUAL
          //==============================================

          for(ii=0;ii<nlbf;ii++)
          {
             bb1 = dN_dx[ii];
             bb2 = dN_dy[ii];

             bc[0][0] = (bb1 * D11[0][0] + bb2 * D11[3][0]);
             bc[0][1] = (bb1 * D11[0][1] + bb2 * D11[3][1]);
             bc[0][2] = (bb1 * D11[0][3] + bb2 * D11[3][3]);

             bc[1][0] = (bb2 * D11[1][0] + bb1 * D11[3][0]);
             bc[1][1] = (bb2 * D11[1][1] + bb1 * D11[3][1]);
             bc[1][2] = (bb2 * D11[1][3] + bb1 * D11[3][3]);

             fact1 = (bb1 * stre[0] + bb2 * stre[3]) ;
             fact2 = (bb1 * stre[3] + bb2 * stre[1]) ;

             TI   = 2*ii;
             TIp1 = TI+1;

             Flocal(TI)   -= fact1 ;
             Flocal(TIp1) -= fact2 ;

             for(jj=0;jj<nlbf;jj++)
             {
                bb1 = dN_dx[jj];
                bb2 = dN_dy[jj];

                TJ   = 2*jj;
                TJp1 = TJ+1;

                fact = fact1 * bb1 + fact2 * bb2;

                Klocal(TI,TJ)     +=  (bc[0][0] * bb1 + bc[0][2] * bb2 + fact) ;
                Klocal(TI,TJp1)   +=  (bc[0][1] * bb2 + bc[0][2] * bb1) ;
                Klocal(TIp1,TJ)   +=  (bc[1][0] * bb1 + bc[1][2] * bb2) ;
                Klocal(TIp1,TJp1) +=  (bc[1][1] * bb2 + bc[1][2] * bb1 + fact) ;
             }
             // compute Kup and Kpp matrices

             bb1 *= dvol;
             bb2 *= dvol;

             for(jj=0;jj<sizep;jj++)
             {
               Kup(TI, jj)   += ( bb1 * Nbar[jj] );
               Kup(TIp1, jj) += ( bb2 * Nbar[jj] );
             }
          }

          fact1 = (detF - 1.0 - pres/BULK) * dvol0;
          fact2 = dvol0/BULK;

          for(ii=0;ii<sizep;ii++)
          {
            //Flocal2(ii) -= Nbar[ii] * fact1;

            fact = Nbar[ii] * fact2;

            for(jj=0;jj<sizep;jj++)
               Kpp(ii,jj) -= ( fact * Nbar[jj] );
          }
    }//gp1
  }//gp2

  return 0;
}


int LagrangeElem3DStructSolidMixed::calcInternalForces()
{
  return 0;
}



void LagrangeElem3DStructSolidMixed::discreteContourplot(int vartype, int varindex, int index, int nCol, double umin, double umax)
{
  return;
}


void LagrangeElem3DStructSolidMixed::projectToKnots(bool extrapolateFlag, int vartype, int varindex, int index)
{
  return;
}


void LagrangeElem3DStructSolidMixed::projectStress(int varindex, double* outval)
{
  return;
}



void LagrangeElem3DStructSolidMixed::projectStrain(int vartype, int varindex, double* outval)
{
  return;
}



void LagrangeElem3DStructSolidMixed::projectIntVar(int index, double* outval)
{
  return;
}


int LagrangeElem3DStructSolidMixed::calcOutput(double u1, double v1)
{
  return 0;
}



void LagrangeElem3DStructSolidMixed::toPostprocess(int vartype, int varindex, int type, SparseMatrixXd&  coeffMat, VectorXd& rhsVec)
{
  return;
}




