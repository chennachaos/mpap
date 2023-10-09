
#include "LagrangeElem2DBbarFbar.h"

#include "Debug.h"
#include "MpapTime.h"
#include "ComputerTime.h"
#include "GeomDataLagrange.h"
#include "SolutionData.h"
#include "FunctionsMaterial.h"
#include "QuadratureUtil.h"
#include "TimeFunction.h"


using namespace std;

extern ComputerTime       computerTime;
extern MpapTime mpapTime;
extern List<TimeFunction> timeFunction;



LagrangeElem2DBbarFbar::LagrangeElem2DBbarFbar()
{
  if (debug) cout << " constructor LagrangeElem2DBbarFbar\n\n";

  ndim   = 2;
  degree = 1;
  npElem = 4;
  nlbf   = 4;
  ndof   = 2;
  nsize  = npElem*ndof;
}


LagrangeElem2DBbarFbar::~LagrangeElem2DBbarFbar()
{
  if (debug) cout << " destructor LagrangeElem2DBbarFbar\n\n";
}


void LagrangeElem2DBbarFbar::prepareElemData()
{
  LagrangeElement::prepareElemData();

  return;
}


void LagrangeElem2DBbarFbar::prepareElemData2()
{
  return;
}


int LagrangeElem2DBbarFbar::calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal, bool firstIter)
{
  if(finite)
    calcStiffnessAndResidualFS(Klocal, Flocal);
  else
    calcStiffnessAndResidualSS(Klocal, Flocal);

  return 0;
}




int LagrangeElem2DBbarFbar::calcStiffnessAndResidualSS(MatrixXd& Klocal, VectorXd& Flocal)
{
    double F[4], Fc[4], detFc, detF=0.0, F33, fact, fact1, fact2, dvol, dvol0;
    double Jac, bb1, bb2, bb3, cc1, cc2, cc3;

    VectorXd  N(nlbf), dN_dx(nlbf), dN_dy(nlbf);
    VectorXd  dispC(nsize), velC(nsize), accC(nsize);
    MatrixXd  Mlocal(nsize, nsize);

    double  stre[4], cc[4][4], bc[2][4], bforce[2];
    double  param[2], Bbar[4][2], r1d3 = 1.0/3.0;

    VectorXd  Wmat1(nlbf), Wmat2(nlbf);

    int   err,  isw,  count,  count1, index, ll = 0, ii, jj, gp, TI, TIp1, TJ, TJp1;
    int   ind1, ind2, kk;

    double af = SolnData->td(2);
    double d1 = SolnData->td(5);
    double aa = SolnData->td(10);
    double dt = mpapTime.dt;

    double xNode[4], yNode[4], xx, yy;

    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    for(ii=0;ii<npElem;ii++)
    {
      ind1 = ndof*ii;
      ind2 = nodeNums[ii]*ndof;

      for(kk=0;kk<ndof;kk++)
      {
        dispC(ind1+kk)  =  SolnData->var1Cur[ind2+kk];
        velC(ind1+kk)   =  SolnData->var1DotCur[ind2+kk];
        accC(ind1+kk)   =  SolnData->var1DotDotCur[ind2+kk];
      }
    }

    if(Klocal.rows() != nsize)
    {
      Klocal.resize(nsize, nsize);
      Flocal.resize(nsize);
    }
    Klocal.setZero();
    Flocal.setZero();
    Mlocal.setZero();

    //MatrixXd  Idev(4,4);
    //double  r1d3 = 1.0/3.0, r2d3 = 2.0*r1d3;
    //Idev[0][0] =  r2d3;     Idev[0][1] = -r1d3;  Idev[0][2] = -r1d3;    Idev[0][3] = 0.0;
    //Idev[1][0] = -r1d3;     Idev[1][1] =  r2d3;  Idev[1][2] = -r1d3;    Idev[1][3] = 0.0;
    //Idev[2][0] = -r1d3;     Idev[2][1] = -r1d3;  Idev[2][2] =  r2d3;    Idev[2][3] = 0.0;
    //Idev[3][0] =  0.0;      Idev[3][1] =  0.0;   Idev[3][2] =  0.0;     Idev[3][3] = 1.0;

    //  compute determinant tr(F) in element centre 

    vector<double>  gausspoints1, gausspoints2, gaussweights;

    getGaussPointsQuad(1, gausspoints1, gausspoints2, gaussweights);

    param[0] = gausspoints1[0];
    param[1] = gausspoints2[0];

    GeomData->computeBasisFunctions2D(0, 2, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

    Fc[0] = computeValueCur(0, dN_dx) + 1.0;
    Fc[2] = computeValueCur(0, dN_dy);
    Fc[1] = computeValueCur(1, dN_dx);
    Fc[3] = computeValueCur(1, dN_dy) + 1.0;

    detFc = Fc[0]*Fc[3] - Fc[1]*Fc[2];

    for(jj=0;jj<nlbf;jj++)
    {
      Wmat1(jj) = dN_dx[jj] ;
      Wmat2(jj) = dN_dy[jj] ;
    }

    double rho0 = elmDat[5] ;
    double rho  = rho0;
    bforce[0]   = elmDat[6]*timeFunction[0].prop ;
    bforce[1]   = elmDat[7]*timeFunction[0].prop ;

    getGaussPointsQuad(nGP, gausspoints1, gausspoints2, gaussweights);

    count = 1;   ll = 0;   err = 0;   isw = 3;
    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(0, 2, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

        fact  = gaussweights[gp] * thick;
        dvol0 = Jac * fact;
        dvol  = dvol0;

        //GeomData->computeDeformationGradient(1, nodeNums, &dN_dx(0), &dN_dy(0), F, detF);

        F[0] = computeValueCur(0, dN_dx) + 1.0;
        F[2] = computeValueCur(0, dN_dy);
        F[1] = computeValueCur(1, dN_dx);
        F[3] = computeValueCur(1, dN_dy) + 1.0;

        detF = F[0]*F[3] - F[1]*F[2];

        //cout << dvol << '\t' << F[0] << '\t' << F[1] << '\t' << F[2] << '\t' << F[3] << endl;

        if(detF < 0.0)   return 1;

        // ADJUST F33 fOR 2D PROBLEMS BASED ON THE ASSUMPTIONS OF PLANE STRESS/PLANE STRAIN/AXISYMMETRIC

        if(sss == 1)  // plane stress
        {
          F33 = 3.0 - F[0] - F[3];
        }
        else if(sss == 2)    // plane strain
          F33 = 1.0;

        fact = (Fc[0] + Fc[3] - F[0] - F[3])/3.0;

        F[0] += fact ;
        F[3] += fact ;
        F33 = 1.0 + fact ;


        matlib2d_(matDat, F, &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);
        count++;
        ll += nivGP;

        if(err !=0)    return 1;

        for(ii=0;ii<4;ii++)
        {
          stre[ii] *= dvol;
          for(jj=0;jj<4;jj++)
            cc[ii][jj] *= dvol;

          //cout << " ii = " << ii << '\t' << stre[ii] << endl;
        }

        //==============================================
        // CALCULATE TANGENT STIFFNESS
        //==============================================

        for(ii=0;ii<nlbf;ii++)
        {
          bb1 = dN_dx[ii];
          bb2 = dN_dy[ii];
          bb3 = N[ii]*dvol0*rho;

          fact1 = (Wmat1(ii) - bb1)/3.0;
          fact2 = (Wmat2(ii) - bb2)/3.0;

          Bbar[0][0] = bb1+fact1;
          Bbar[1][0] = fact1;
          Bbar[2][0] = fact1;
          Bbar[3][0] = bb2;

          Bbar[0][1] = fact2;
          Bbar[1][1] = bb2+fact2;
          Bbar[2][1] = fact2;
          Bbar[3][1] = bb1;

          //Bbar[0][0] = bb1;
          //Bbar[1][0] = -r1d3*bb1;
          //Bbar[2][0] = -r1d3*bb1;
          //Bbar[3][0] = bb2;

          //Bbar[0][1] = -r1d3*bb2;
          //Bbar[1][1] = bb2;
          //Bbar[2][1] = -r1d3*bb2;
          //Bbar[3][1] = bb1;

          bc[0][0] = (Bbar[0][0] * cc[0][0] + Bbar[1][0] * cc[1][0] + Bbar[2][0] * cc[2][0] + Bbar[3][0] * cc[3][0]);
          bc[0][1] = (Bbar[0][0] * cc[0][1] + Bbar[1][0] * cc[1][1] + Bbar[2][0] * cc[2][1] + Bbar[3][0] * cc[3][1]);
          bc[0][2] = (Bbar[0][0] * cc[0][2] + Bbar[1][0] * cc[1][2] + Bbar[2][0] * cc[2][2] + Bbar[3][0] * cc[3][2]);
          bc[0][3] = (Bbar[0][0] * cc[0][3] + Bbar[1][0] * cc[1][3] + Bbar[2][0] * cc[2][3] + Bbar[3][0] * cc[3][3]);

          bc[1][0] = (Bbar[0][1] * cc[0][0] + Bbar[1][1] * cc[1][0] + Bbar[2][1] * cc[2][0] + Bbar[3][1] * cc[3][0]);
          bc[1][1] = (Bbar[0][1] * cc[0][1] + Bbar[1][1] * cc[1][1] + Bbar[2][1] * cc[2][1] + Bbar[3][1] * cc[3][1]);
          bc[1][2] = (Bbar[0][1] * cc[0][2] + Bbar[1][1] * cc[1][2] + Bbar[2][1] * cc[2][2] + Bbar[3][1] * cc[3][2]);
          bc[1][3] = (Bbar[0][1] * cc[0][3] + Bbar[1][1] * cc[1][3] + Bbar[2][1] * cc[2][3] + Bbar[3][1] * cc[3][3]);

          TI   = 2*ii;
          TIp1 = TI+1;

          Flocal(TI)   += bb3*bforce[0];
          Flocal(TIp1) += bb3*bforce[1];

          Flocal(TI)   -= (Bbar[0][0]*stre[0] + Bbar[1][0]*stre[1] + Bbar[2][0]*stre[2] + Bbar[3][0]*stre[3]) ;
          Flocal(TIp1) -= (Bbar[0][1]*stre[0] + Bbar[1][1]*stre[1] + Bbar[2][1]*stre[2] + Bbar[3][1]*stre[3]) ;

          for(jj=0; jj<nlbf; jj++)
          {
            cc1 = dN_dx[jj];
            cc2 = dN_dy[jj];

            fact1 = (Wmat1(jj) - cc1)/3.0;
            fact2 = (Wmat2(jj) - cc2)/3.0;

            Bbar[0][0] = cc1+fact1;
            Bbar[1][0] = fact1;
            Bbar[2][0] = fact1;
            Bbar[3][0] = cc2;

            Bbar[0][1] = fact2;
            Bbar[1][1] = cc2+fact2;
            Bbar[2][1] = fact2;
            Bbar[3][1] = cc1;

            TJ   = 2*jj;
            TJp1 = TJ+1;

            Klocal(TI,TJ)     += af*(bc[0][0] * Bbar[0][0] + bc[0][1] * Bbar[1][0] + bc[0][2] * Bbar[2][0] + bc[0][3] * Bbar[3][0]) ;
            Klocal(TI,TJp1)   += af*(bc[0][0] * Bbar[0][1] + bc[0][1] * Bbar[1][1] + bc[0][2] * Bbar[2][1] + bc[0][3] * Bbar[3][1]) ;
            Klocal(TIp1,TJ)   += af*(bc[1][0] * Bbar[0][0] + bc[1][1] * Bbar[1][0] + bc[1][2] * Bbar[2][0] + bc[1][3] * Bbar[3][0]) ;
            Klocal(TIp1,TJp1) += af*(bc[1][0] * Bbar[0][1] + bc[1][1] * Bbar[1][1] + bc[1][2] * Bbar[2][1] + bc[1][3] * Bbar[3][1]) ;

            Mlocal(TI,  TJ)   +=  bb3*N[jj];
            Mlocal(TIp1,TJp1) +=  bb3*N[jj];
          }
        }
    }//gp

    //printMatrix(Klocal);
    //if(elenum == 10)
      //printMatrix(Klocal);  printf("\n\n\n");  printMatrix(Mlocal);  printf("\n\n\n"); printVector(Flocal);

    Klocal +=  d1*Mlocal;
    Flocal -=  Mlocal*accC;

    //Klocal /= aa;
    //cout << '\t' << " Total Volume for element # " << elenum << " is = " << totvol << endl; cout << endl;

    return 0;
}




int LagrangeElem2DBbarFbar::calcStiffnessAndResidualFS(MatrixXd& Klocal, VectorXd& Flocal)
{
    double  F[4], Fc[4], detF=0.0, F33, fact, fact1, fact2, dvol, dvol0, Jac;
    double  bb1, bb2, bb3, cc1, cc2, cc3;

    VectorXd  N(nlbf), dN_dx(nlbf), dN_dy(nlbf);
    VectorXd  dispC(nsize), velC(nsize), accC(nsize);
    MatrixXd  Mlocal(nsize, nsize);

    double  stre[4], cc[4][4], cc0[4][4], bc[2][4], param[2], cch[4], bforce[2];
    double  detFc, trFc;
    VectorXd  Wmat1(nlbf), Wmat2(nlbf);

    int  ind1, ind2, kk, err, isw, count, count1, index, ll = 0, ii, jj, gp, TI, TIp1, TJ, TJp1;

    double rho0 = elmDat[5] ;
    double rho  = rho0;
    bforce[0]   = rho0*elmDat[6] ;
    bforce[1]   = rho0*elmDat[7] ;

    double af = SolnData->td(2);
    double d1 = SolnData->td(5);
    double aa = SolnData->td(10);
    double dt = mpapTime.dt;

    double xNode[4], yNode[4], xx, yy;

    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    for(ii=0;ii<npElem;ii++)
    {
      ind1 = ndof*ii;
      ind2 = nodeNums[ii]*ndof;

      for(kk=0;kk<ndof;kk++)
      {
        dispC(ind1+kk)  =  SolnData->var1Cur[ind2+kk];
        velC(ind1+kk)   =  SolnData->var1DotCur[ind2+kk];
        accC(ind1+kk)   =  SolnData->var1DotDotCur[ind2+kk];
      }
    }

    // resize local matrices and initialise them to zero
    if(Klocal.rows() != nsize)
    {
      Klocal.resize(nsize, nsize);
      Flocal.resize(nsize);
    }
    Klocal.setZero();
    Flocal.setZero();
    Mlocal.setZero();

    //  compute determinant detF (tr(F)) in element center
    //
    vector<double>  gausspoints1, gausspoints2, gaussweights;

    getGaussPointsQuad(1, gausspoints1, gausspoints2, gaussweights);

    param[0] = gausspoints1[0];
    param[1] = gausspoints2[0];

    GeomData->computeBasisFunctions2D(0, 2, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

    Fc[0] = computeValueCur(0, dN_dx) + 1.0;
    Fc[2] = computeValueCur(0, dN_dy);
    Fc[1] = computeValueCur(1, dN_dx);
    Fc[3] = computeValueCur(1, dN_dy) + 1.0;

    detFc = Fc[0]*Fc[3] - Fc[1]*Fc[2];
    trFc  = Fc[0] + Fc[3];

    //cout << " trFc = " << trFc << '\t' << detFc << '\t' << Jac << endl;

    if(finite)
      GeomData->computeBasisFunctions2D(1, 2, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

    for(jj=0;jj<nlbf;jj++)
    {
      Wmat1(jj) = dN_dx[jj] ;
      Wmat2(jj) = dN_dy[jj] ;
    }


    getGaussPointsQuad(nGP, gausspoints1, gausspoints2, gaussweights);

    count = 1;   ll = 0;   err = 0;   isw = 3;
    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(0, 2, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

        //for(ii=0; ii<nlbf; ii++)
          //cout << ii << '\t' << dN_dx[ii] << '\t' << dN_dy[ii] << endl;

        fact  = gaussweights[gp] * thick;
        dvol0 = Jac * fact;
        dvol  = dvol0;

        //GeomData->computeDeformationGradient(1, nodeNums, &dN_dx(0), &dN_dy(0), F, detF);

        F[0] = computeValueCur(0, dN_dx) + 1.0;
        F[2] = computeValueCur(0, dN_dy);
        F[1] = computeValueCur(1, dN_dx);
        F[3] = computeValueCur(1, dN_dy) + 1.0;

        detF = F[0]*F[3] - F[1]*F[2];

        //printf("F... \t%20.18f\t%20.18f\t%20.18f\t%20.18f\t%20.18f \n", F[0], F[1], F[2], F[3], detF);
        //cout << " finite = " << finite << endl;

        if(detF < 0.0)
          return 1;

        if(finite)
        {
          GeomData->computeBasisFunctions2D(1, 2, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);
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

        //    replace detF (tr(F)) with detF (tr(F)) from element centre

        if(finite)
        {
          fact = sqrt(detFc/detF);

          for(ii=0; ii<4; ii++)
            F[ii] = F[ii] * fact ;
        }
        else
        {
            fact = 0.5 * (trFc - F[0] - F[3]) ;
            F[0] = F[0] + fact ;
            F[3] = F[3] + fact ;
        }

        matlib2d_(matDat, F, &F33, stre, cc0[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);
        count++;
        ll += nivGP;

        //printf(" stress... \t%14.12f\t%14.12f\t%14.12f\t%14.12f\n\n", stre[0], stre[1], stre[2], stre[3]);
        //printf(" dvol = %14.12f \n", dvol);

        if(err !=0)
          return 1;

        for(ii=0;ii<4;ii++)
        {
          for(jj=0;jj<4;jj++)
          {
            cc[ii][jj] = cc0[ii][jj] ;
          }
        }

        /*
        cc[0][2] = cc1[0][3];
        cc[1][2] = cc1[1][3];

        cc[2][0] = cc1[3][0];
        cc[2][1] = cc1[3][1];

        cc[2][2] = cc1[3][3];
        cc[2][3] = cc1[3][3];
        cc[3][2] = cc1[3][3];
        */

        if(finite)
          dvol *= F33;

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

        // correct tangent cc for Fbar forumation
        //cch[0]  = 0.5 * (cc[0][0] + cc[0][1] - stre[0]);
        //cch[1]  = 0.5 * (cc[1][0] + cc[1][1] - stre[1]);
        //cch[2]  = 0.5 * (cc[3][0] + cc[3][1] - stre[3]);

        cch[0]  = 0.5 * (cc[0][0] + cc[0][1] );
        cch[1]  = 0.5 * (cc[1][0] + cc[1][1] );
        cch[2]  = 0.5 * (cc[3][0] + cc[3][1] );

        cc[0][0] = cc[0][0] - cch[0];
        cc[1][0] = cc[1][0] - cch[1];
        cc[3][0] = cc[3][0] - cch[2];

        cc[0][1] = cc[0][1] - cch[0];
        cc[1][1] = cc[1][1] - cch[1];
        cc[3][1] = cc[3][1] - cch[2];

        for(ii=0;ii<nlbf;ii++)
        {
          bb1 = dN_dx[ii];
          bb2 = dN_dy[ii];
          bb3 = N[ii]*dvol0*rho;

          fact1 = bb1 * cch[0] + bb2 * cch[2] ;
          fact2 = bb1 * cch[2] + bb2 * cch[1] ;

          bc[0][0] = bb1 * cc[0][0] + bb2 * cc[3][0];
          bc[0][1] = bb1 * cc[0][1] + bb2 * cc[3][1];
          bc[0][2] = bb1 * cc[0][3] + bb2 * cc[3][3];

          bc[1][0] = bb2 * cc[1][0] + bb1 * cc[3][0];
          bc[1][1] = bb2 * cc[1][1] + bb1 * cc[3][1];
          bc[1][2] = bb2 * cc[1][3] + bb1 * cc[3][3];

          TI   = 2*ii;
          TIp1 = TI+1;

          Flocal(TI)   += bb3*bforce[0];
          Flocal(TIp1) += bb3*bforce[1];

          Flocal(TI)   -= (bb1*stre[0] + bb2*stre[3]) ;
          Flocal(TIp1) -= (bb1*stre[3] + bb2*stre[1]) ;

          for(jj=0; jj<nlbf; jj++)
          {
              cc1 = dN_dx[jj];
              cc2 = dN_dy[jj];

              TJ   = 2*jj;
              TJp1 = TJ+1;

              Klocal(TI,   TJ)   +=  af*(bc[0][0] * cc1 + bc[0][2] * cc2) ;
              Klocal(TI,   TJp1) +=  af*(bc[0][1] * cc2 + bc[0][2] * cc1) ;
              Klocal(TIp1, TJ)   +=  af*(bc[1][0] * cc1 + bc[1][2] * cc2) ;
              Klocal(TIp1, TJp1) +=  af*(bc[1][1] * cc2 + bc[1][2] * cc1) ;

              //Klocal(TI,TJ)     +=  ( fact1 * (Wmat1(jj)-bb1) );
              //Klocal(TI,TJp1)   +=  ( fact1 * (Wmat2(jj)-bb2) );
              //Klocal(TIp1,TJ)   +=  ( fact2 * (Wmat1(jj)-bb1) );
              //Klocal(TIp1,TJp1) +=  ( fact2 * (Wmat2(jj)-bb2) );

              Klocal(TI,   TJ)   +=  af*( fact1 * Wmat1(jj) );
              Klocal(TI,   TJp1) +=  af*( fact1 * Wmat2(jj) );
              Klocal(TIp1, TJ)   +=  af*( fact2 * Wmat1(jj) );
              Klocal(TIp1, TJp1) +=  af*( fact2 * Wmat2(jj) );

              Mlocal(TI,   TJ)   +=  bb3*N[jj];
              Mlocal(TIp1, TJp1) +=  bb3*N[jj];
          }
        }

        //   part 2. -- geometrical matrix  (if geometry nonlinear)

        if(finite)
        {
          for(ii=0; ii<nlbf; ii++)
          {
              fact1 = dN_dx[ii]*stre[0] + dN_dy[ii]*stre[3] ;
              fact2 = dN_dx[ii]*stre[3] + dN_dy[ii]*stre[1] ;

              TI   = 2*ii;
              TIp1 = TI+1;

              for(jj=0; jj<nlbf; jj++)
              {
                TJ   = 2*jj;
                TJp1 = TJ+1;

                fact = af*(fact1*dN_dx[jj] + fact2*dN_dy[jj]);

                Klocal(TI,   TJ)   += fact ;
                Klocal(TIp1, TJp1) += fact ;
              }
          }
        } // if(finite)
    } //gp

    Klocal +=  d1*Mlocal;
    Flocal -=  Mlocal*accC;

    //Klocal /= aa;

    /*
    int vec[]={0,1,2,3,6,7,4,5};

    MatrixXd  K2(nsize, nsize);
    for(ii=0; ii<nsize; ii++)
    {
      for(jj=0; jj<nsize; jj++)
        K2(ii,jj) = Klocal(vec[ii], vec[jj]);
    }

    printMatrix(Klocal);
    printf("\n\n\n");
    printMatrix(K2);
    //printVector(Flocal);
    */

    return 0;
}



int LagrangeElem2DBbarFbar::calcLoadVector()
{
  return 0;
}



int LagrangeElem2DBbarFbar::calcInternalForces()
{

  return 0;
}



void LagrangeElem2DBbarFbar::discreteContourplot(int vartype, int varindex, int index, int nCol, double umin, double umax)
{
  return;
}


void LagrangeElem2DBbarFbar::projectToKnots(bool extrapolateFlag, int vartype, int varindex, int index)
{
  return;
}


void LagrangeElem2DBbarFbar::projectStress(int varindex, double* outval)
{
  return;
}



void LagrangeElem2DBbarFbar::projectStrain(int vartype, int varindex, double* outval)
{
  return;
}



void LagrangeElem2DBbarFbar::projectIntVar(int index, double* outval)
{
  return;
}


int LagrangeElem2DBbarFbar::calcOutput(double u1, double v1)
{
  return 0;
}



void LagrangeElem2DBbarFbar::toPostprocess(int vartype, int varindex, int type, SparseMatrixXd&  coeffMat, VectorXd& rhsVec)
{

  return;
}


void  LagrangeElem2DBbarFbar::computeEnergy(int index1, int index2, VectorXd& energy)
{
    // compute total energy for structural dynamic problems
    ///////////////////////////////////////////////////////////

    elmDat = &(SolnData->ElemProp[elmType].data[0]);
    //matDat = &(SolidSolnData->MatlProp[matType].data[0]);

  double  F[4], Fc[4], detF, detFc, trFc, F33, fact, dvol, dvol0, Jac, dt, JacMult;
  double  posi[2], disp[2], velo[2];
  double  stress[4], strain[4], cc[4][4], param[2], bforce[2];

  VectorXd  N(nlbf), dN_dx(nlbf), dN_dy(nlbf);

  int   err,  isw,  count,  count1, index, ll = 0, ii, jj, gp1, gp2, TI, TIp1, TJ, TJp1;
  int   ind1, ind2, kk;
  
  double  rho  = elmDat[5] ;

  vector<double>  gausspoints1, gaussweights1;

  //  compute determinant detF (tr(F)) in element centre 

  getGaussPoints1D(1, gausspoints1, gaussweights1);

  param[0] = GeomData->gausspoints1[0];
  param[1] = GeomData->gausspoints1[0];

  GeomData->computeBasisFunctions2D(0, 2, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

  Fc[0] = computeValue(0, dN_dx) + 1.0;
  Fc[2] = computeValue(0, dN_dy);
  Fc[1] = computeValue(1, dN_dx);
  Fc[3] = computeValue(1, dN_dy) + 1.0;

  detFc = Fc[0]*Fc[3] - Fc[1]*Fc[2];
  trFc  = Fc[0] + Fc[3];



  count = 1;   ll = 0;   err = 0;   isw = 3;
  dt = mpapTime.dt;

  //cout << " finite = " << finite << endl;

  int nGP = (int) elmDat[0] ;

  getGaussPoints1D(nGP, gausspoints1, gaussweights1);

  energy.setZero();

  for(gp2=0;gp2<nGP;gp2++)
  {
    JacMult = gaussweights1[gp2] * thick;

    param[1] = gausspoints1[gp2];

  for(gp1=0;gp1<nGP;gp1++)
  {
        param[0] = gausspoints1[gp1];

        GeomData->computeBasisFunctions2D(0, 2, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

        fact = gaussweights1[gp1] * JacMult;

        dvol0 = Jac * fact;
        dvol = dvol0;

        //if(axsy)
          //dvol *= 2.0*PI*yy;

        //GeomData->computeDeformationGradient(1, nodeNums, &dN_dx(0), &dN_dy(0), F, detF);

        F[0] = computeValue(0, dN_dx) + 1.0;
        F[2] = computeValue(0, dN_dy);
        F[1] = computeValue(1, dN_dx);
        F[3] = computeValue(1, dN_dy) + 1.0;

        detF =  F[0]*F[3] - F[1]*F[2];

        //if(detF < 0.0)   return 1;

        //if(finite)
        //{
          //GeomData->computeBasisFunctions2D(1, 2, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);
          //dvol = Jac * fact;
        //}

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
        
        //printf(" %14.12f \t %14.12f \t %14.12f \t %14.12f \n ", F[0], F[1], F[2], F[3]);

        if(finite)
        {
          fact = sqrt(detFc/detF);

          for(ii=0; ii<4; ii++)
            F[ii] = F[ii] * fact ;
        }
        else
        {
            fact = 0.5 * (trFc - F[0] - F[3]) ;
            F[0] = F[0] + fact ;
            F[3] = F[3] + fact ;
        }

        matlib2d_(matDat, F, &F33, stress, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);
        count++;
        ll += nivGP;

        strain[0] = F[0] - 1.0;
        strain[1] = F[3] - 1.0;
        strain[2] = F33  - 1.0;
        strain[3] = 0.5 * (F[1] + F[2]);

        //printf(" %14.12f \t %14.12f \t %14.12f \t %14.12f \n ", strain[0], strain[1], strain[2], strain[3]);
        //printf(" %14.12f \t %14.12f \t %14.12f \t %14.12f \n ", stress[0], stress[1], stress[2], stress[3]);
        //printf("  dvol =  %14.12f \n\n ", dvol);

        velo[0] = computeValueDot(0, N);
        velo[1] = computeValueDot(1, N);

        // kinetic energy
        energy[0] += (0.5*dvol0*rho*(velo[0]*velo[0]+velo[1]*velo[1]));
        energy[1] += (0.5*dvol*(strain[0]*stress[0]+strain[1]*stress[1]+strain[2]*stress[2]+strain[3]*stress[3]));

  }//gp1
  }//gp2

  //printf("\n\n");

  return;
}





