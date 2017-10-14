
#include "LagrangeElem3DStructSolidTet4NodeStab.h"

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



LagrangeElem3DStructSolidTet4NodeStab::LagrangeElem3DStructSolidTet4NodeStab()
{
  if (debug) cout << " constructor LagrangeElem3DStructSolidTet4NodeStab\n\n";

  ndim   = 3;
  degree = 1;
  npElem = 4;
  nlbf   = 4;
  ndof   = 4;
  nsize  = nlbf*ndof;

}

LagrangeElem3DStructSolidTet4NodeStab::~LagrangeElem3DStructSolidTet4NodeStab()
{
  if (debug) cout << " destructor LagrangeElem3DStructSolidTet4NodeStab\n\n";
}


void LagrangeElem3DStructSolidTet4NodeStab::prepareElemData()
{
  LagrangeElement::prepareElemData();

  int ii, jj;

  ndim   = 3;
  degree = 1;
  npElem = 4;
  nlbf   = 4;
  ndof   = 4;
  nsize  = nlbf*ndof;


  // set the element property variables

  elmDat = &(SolnData->ElemProp[elmType].data[0]);
  matDat = &(SolnData->MatlProp[matType].data[0]);

  finiteInt  = (int) elmDat[2] ;

  matId      = SolnData->MatlProp[matType].id + 1;
  finite     = (finiteInt >= 1) ;

  //Klocal.resize(nsize, nsize);
  //Flocal.resize(nsize);
  
  return;
}



void LagrangeElem3DStructSolidTet4NodeStab::prepareElemData2()
{
  cout << " LagrangeElem3DStructSolidTet4NodeStab::prepareElemData2() " << endl;
  return;
}

int LagrangeElem3DStructSolidTet4NodeStab::calcLoadVector()
{
  return 0;
}


int LagrangeElem3DStructSolidTet4NodeStab::calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal, bool firstIter)
{
//   char fct[] = "LagrangeElem3DStructSolidTet4NodeStab::calcStiffnessAndResidual";
//   computerTime.go(fct);

  double  F[9], detF=0.0, F33, fact, fact1, fact2, dvol, dvol0, tau, r2d3=2.0/3.0;
  double  Jac, JacTemp, dt, totvol, volstr, pres, pbar, dp[3];
  double  bb1, bb2, bb3, bb4, bb5, cc1, cc2, cc3, cc4, cc5;

  VectorXd  N(nlbf), dN_dx(nlbf), dN_dy(nlbf), dN_dz(nlbf);
  MatrixXd  Mlocal(nsize, nsize), Dj(3,4);
  MatrixXd  forVolume(4,4);

  double  stre[6], strdev[6], cc[6][6], cctmp[6][6], bc[3][6], param[3], force[3];
  double  bforce[3], sig[3], acceCur[3], rStab[3], Idev[6][6];

  Idev3D(Idev);


  int   err,  isw,  count,  count1, index, ll = 0, ii, jj, kk, mm, gp;
  int   TI, TIp1, TIp2, TIp3, TJ, TJp1, TJp2, TJp3;
  int   ind1, ind2;
  
  double  BULK = matDat[0] ;
  double  mu   = matDat[1] ;
  double  finiteFact = finite ? 1.0 : 0.0 ;
  double  eps  = 1.0/BULK;

  double rho0 = elmDat[5] ;
  bforce[0]  = 0.0;
  bforce[1]  = 0.0;
  bforce[2]  = 0.0;

  //bforce[0]  = elmDat[6]*timeFunction[0].prop ;
  //bforce[1]  = elmDat[7]*timeFunction[0].prop ;
  //bforce[2]  = elmDat[7]*timeFunction[0].prop ;

  double  af = SolnData->td(2);
  double  acceFact1 = SolnData->td(5);
  double  acceFact2 = acceFact1;

  //printf(" %14.12f \t %14.12f \t %14.12f \t %14.12f \n ", BULK, mu, bforce[0], bforce[1]);
  //printf(" %14.12f \t %14.12f \t %14.12f \t %14.12f \n ", rho, af, d1, aa);
  //printf(" %5d \t %5d \t %5d \t %5d \n ", nodeNums[0], nodeNums[1], nodeNums[2], nodeNums[3]);

  double xNode[4], yNode[4], zNode[4], xx, yy, zz;

  for(ii=0;ii<npElem;ii++)
  {
    xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
    yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    zNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][2];
  }

  //forVolume(0,0) = 1.0; forVolume(0,1) = xNode[0]; forVolume(0,2) = yNode[0]; forVolume(0,3) = zNode[0];
  //forVolume(1,0) = 1.0; forVolume(1,1) = xNode[1]; forVolume(1,2) = yNode[1]; forVolume(1,3) = zNode[1];
  //forVolume(2,0) = 1.0; forVolume(2,1) = xNode[2]; forVolume(2,2) = yNode[2]; forVolume(2,3) = zNode[2];
  //forVolume(3,0) = 1.0; forVolume(3,1) = xNode[3]; forVolume(3,2) = yNode[3]; forVolume(3,3) = zNode[3];

  forVolume(0,0) = 1.0;      forVolume(0,1) = 1.0;      forVolume(0,2) = 1.0;      forVolume(0,3) = 1.0;
  forVolume(1,0) = xNode[0]; forVolume(1,1) = xNode[1]; forVolume(1,2) = xNode[2]; forVolume(1,3) = xNode[3];
  forVolume(2,0) = yNode[0]; forVolume(2,1) = yNode[1]; forVolume(2,2) = yNode[2]; forVolume(2,3) = yNode[3];
  forVolume(3,0) = zNode[0]; forVolume(3,1) = zNode[1]; forVolume(3,2) = zNode[2]; forVolume(3,3) = zNode[3];

  totvol = forVolume.determinant()/6.0;
  //cout << " totvol = " << totvol << endl;

  //double  h2  = totvol*4.0/PI;
  //double h = pow(3.0*totvol/PI/4.0, 1.0/3.0);
  double h = 1.0/6.0;
  h = sqrt(3.0*h*h);
  
  double  alpha = 1.0/12.0;

  //tau = elmDat[8]*alpha*h*h/sqrt(mu);
  tau = elmDat[8]*alpha*h*h/mu/dt;

  //printf(" %14.12f \t %14.12f \t %14.12f \t %14.12f \n ", xNode[0], xNode[1], xNode[2], xNode[3]);
  //printf(" %14.12f \t %14.12f \t %14.12f \t %14.12f \n\n\n ", yNode[0], yNode[1], yNode[2], yNode[3]);

    if(Klocal.rows() != nsize)
    {
      Klocal.resize(nsize, nsize);
      Flocal.resize(nsize);
    }
    Klocal.setZero();
    Flocal.setZero();
  Mlocal.setZero();


  count = 1;   ll = 0;   err = 0;   isw = 3;
  dt = mpapTime.dt;

  //cout << " finite = " << finite << endl;

  nGP = (int) elmDat[0] ;

  vector<double>  gausspoints1, gausspoints2, gausspoints3, gaussweights;

  getGaussPointsTet(nGP, gausspoints1, gausspoints2, gausspoints3, gaussweights);


  totvol = 0.0;
  for(gp=0;gp<nGP;gp++)
  {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];
        param[2] = gausspoints3[gp];

        GeomData->computeBasisFunctions3D(0, 1, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);

        JacTemp = gaussweights[gp];

        dvol0 = Jac * JacTemp;
        dvol = dvol0;
        totvol += dvol;

        //cout << " dvol = " << dvol << endl;

        //for(ii=0; ii<nlbf; ii++)
          //cout << N[ii] << '\t' << dN_dx[ii] << '\t' << dN_dy[ii] << endl;

        //GeomData->computeDeformationGradient(1, nodeNums, &dN_dx(0), &dN_dy(0), F, detF);

        F[0] = computeValueCur(0, dN_dx) + 1.0;
        F[3] = computeValueCur(0, dN_dy);
        F[6] = computeValueCur(0, dN_dz);

        F[1] = computeValueCur(1, dN_dx);
        F[4] = computeValueCur(1, dN_dy) + 1.0;
        F[7] = computeValueCur(1, dN_dz);

        F[2] = computeValueCur(2, dN_dx);
        F[5] = computeValueCur(2, dN_dy);
        F[8] = computeValueCur(2, dN_dz) + 1.0;

        acceCur[0] = computeValueDotDotCur(0, N);
        acceCur[1] = computeValueDotDotCur(1, N);
        acceCur[2] = computeValueDotDotCur(2, N);

        pres  = computeValueCur(3, N);

        detF = F[0]*(F[4]*F[8] - F[5]*F[7]) - F[3]*(F[1]*F[8] - F[2]*F[7]) + F[6]*(F[1]*F[5] - F[2]*F[4]);

        if(detF < 0.0)   return 1;

        //printf(" %14.12f \t %14.12f \t %14.12f \t %14.12f \n ", F[0], F[1], F[2], F[3]);

        matlib3d_(matDat, F, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &isw, &err, &count, NULL);
        count++;
        ll += nivGP;

        if(err !=0)    return 1;

        // modify the D matrix
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

        // modify stress vector
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


        if(finite)
        {
          GeomData->computeBasisFunctions3D(1, 1, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);
          dvol = Jac * JacTemp;

          fact  = pbar - pres;
          fact1 = r2d3*pbar - pres;
          fact2 = 2.0*fact;

          for(ii=0;ii<3;ii++)
          {
            for(jj=0;jj<6;jj++)
            {
              cc[ii][jj] -= r2d3*strdev[jj] ;
              cc[jj][ii] -= r2d3*strdev[jj] ;
            }
            for(jj=0;jj<3;jj++)
              cc[ii][jj] -= fact1;

            cc[ii][ii] += fact2;
            kk = ii+3;
            cc[kk][kk] += fact;
          }

          volstr = detF - 1.0 - eps*pres;
        }
        else
        {
          volstr = F[0]+F[4]+F[8]-3.0 - eps*pres;
        }

        //==============================================
        // CALCULATE TANGENT STIFFNESS
        //==============================================

        dp[0] = computeValueCur(3, dN_dx);
        dp[1] = computeValueCur(3, dN_dy);
        dp[2] = computeValueCur(3, dN_dz);

        force[0] = bforce[0]*rho0;
        force[1] = bforce[1]*rho0;
        force[2] = bforce[2]*rho0;

        force[0] -= rho0*acceCur[0];
        force[1] -= rho0*acceCur[1];
        force[2] -= rho0*acceCur[2];

        rStab[0] = - dp[0] - force[0];
        rStab[1] = - dp[1] - force[1];
        rStab[2] = - dp[2] - force[2];

        for(ii=0;ii<nlbf;ii++)
        {
          bb1 = dN_dx[ii]*dvol;
          bb2 = dN_dy[ii]*dvol;
          bb3 = dN_dz[ii]*dvol;
          bb4 = N[ii]*dvol;
          bb5 = N[ii]*dvol0;

          for(kk=0;kk<6;kk++)
          {
            bc[0][kk] = (bb1 * cc[0][kk] + bb2 * cc[3][kk] + bb3 * cc[5][kk]);
            bc[1][kk] = (bb1 * cc[3][kk] + bb2 * cc[1][kk] + bb3 * cc[4][kk]);
            bc[2][kk] = (bb1 * cc[5][kk] + bb2 * cc[4][kk] + bb3 * cc[2][kk]);
          }

          TI   = ndof*ii;
          TIp1 = TI+1;
          TIp2 = TI+2;
          TIp3 = TI+3;

          sig[0] = bb1*stre[0] + bb2*stre[3] + bb3*stre[5] ;
          sig[1] = bb1*stre[3] + bb2*stre[1] + bb3*stre[4] ;
          sig[2] = bb1*stre[5] + bb2*stre[4] + bb3*stre[2] ;

          Flocal(TI)   += (bb5*force[0] - sig[0]) ;
          Flocal(TIp1) += (bb5*force[1] - sig[1]) ;
          Flocal(TIp2) += (bb5*force[2] - sig[2]) ;
          Flocal(TIp3) -= (bb5*volstr)*rho0;

          // PSPG stabilization
          Flocal(TIp3) -= tau*(bb1*rStab[0] + bb2*rStab[1] + bb3*rStab[2]);

          for(jj=0; jj<nlbf; jj++)
          {
              cc1 = dN_dx[jj];
              cc2 = dN_dy[jj];
              cc3 = dN_dz[jj];
              cc4 = N[jj];

              TJ   = ndof*jj;
              TJp1 = TJ+1;
              TJp2 = TJ+2;
              TJp3 = TJ+3;

              acceFact2 = acceFact1*cc4*rho0;

              fact  = bb5*acceFact2;
              fact += af*(sig[0]*cc1 + sig[1]*cc2 + sig[2]*cc3)*finiteFact;

              Klocal(TI,   TJ)    +=  fact;
              Klocal(TIp1, TJp1)  +=  fact;
              Klocal(TIp2, TJp2)  +=  fact;

              Klocal(TI,   TJ)    +=  af*(bc[0][0] * cc1 + bc[0][3] * cc2 + bc[0][5] * cc3) ;
              Klocal(TI,   TJp1)  +=  af*(bc[0][1] * cc2 + bc[0][3] * cc1 + bc[0][4] * cc3) ;
              Klocal(TI,   TJp2)  +=  af*(bc[0][2] * cc3 + bc[0][4] * cc2 + bc[0][5] * cc1) ;
              Klocal(TI,   TJp3)  +=  af*(bb1 * cc4) ;

              Klocal(TIp1, TJ)    +=  af*(bc[1][0] * cc1 + bc[1][3] * cc2 + bc[1][5] * cc3) ;
              Klocal(TIp1, TJp1)  +=  af*(bc[1][1] * cc2 + bc[1][3] * cc1 + bc[1][4] * cc3) ;
              Klocal(TIp1, TJp2)  +=  af*(bc[1][2] * cc3 + bc[1][4] * cc2 + bc[1][5] * cc1) ;
              Klocal(TIp1, TJp3)  +=  af*(bb2 * cc4) ;

              Klocal(TIp2, TJ)    +=  af*(bc[2][0] * cc1 + bc[2][3] * cc2 + bc[2][5] * cc3) ;
              Klocal(TIp2, TJp1)  +=  af*(bc[2][1] * cc2 + bc[2][3] * cc1 + bc[2][4] * cc3) ;
              Klocal(TIp2, TJp2)  +=  af*(bc[2][2] * cc3 + bc[2][4] * cc2 + bc[2][5] * cc1) ;
              Klocal(TIp2, TJp3)  +=  af*(bb3 * cc4) ;

              Klocal(TIp3, TJ)    +=  af*(bb4 * cc1)*rho0 ;
              Klocal(TIp3, TJp1)  +=  af*(bb4 * cc2)*rho0 ;
              Klocal(TIp3, TJp2)  +=  af*(bb4 * cc3)*rho0 ;
              Klocal(TIp3, TJp3)  -=  af*(bb5 * cc4)*eps*rho0;

              // PSPG stabilization

              //fact2 -= ( muTaf*d2N(jj) );
              //fact2 = 0.0;
              fact2 = acceFact2;

              Dj(0,0) = fact2;
              Dj(0,1) = 0.0;
              Dj(0,2) = 0.0;
              Dj(0,3) = -af*dN_dx[jj];

              Dj(1,0) = 0.0;
              Dj(1,1) = fact2;
              Dj(1,2) = 0.0;
              Dj(1,3) = -af*dN_dy[jj];

              Dj(2,0) = 0.0;
              Dj(2,1) = 0.0;
              Dj(2,2) = fact2;
              Dj(2,3) = -af*dN_dz[jj];

              Klocal(TIp3, TJ)   += (bb1*Dj(0,0) + bb2*Dj(1,0) + bb3*Dj(2,0))*tau;
              Klocal(TIp3, TJp1) += (bb1*Dj(0,1) + bb2*Dj(1,1) + bb3*Dj(2,1))*tau;
              Klocal(TIp3, TJp2) += (bb1*Dj(0,2) + bb2*Dj(1,2) + bb3*Dj(2,2))*tau;
              Klocal(TIp3, TJp3) += (bb1*Dj(0,3) + bb2*Dj(1,3) + bb3*Dj(2,3))*tau;
          }
        }
  }//gp
  //cout << " totvol = " << totvol << endl;
  //printMatrix(Klocal);  printf("\n\n\n");  printVector(Flocal);
  //printVector(Flocal);


  return 0;
}






int LagrangeElem3DStructSolidTet4NodeStab::calcInternalForces()
{
  return 0;
}



void LagrangeElem3DStructSolidTet4NodeStab::discreteContourplot(int vartype, int varindex, int index, int nCol, double umin, double umax)
{
  return;
}


void LagrangeElem3DStructSolidTet4NodeStab::projectToKnots(bool extrapolateFlag, int vartype, int varindex, int index)
{
  return;
}


void LagrangeElem3DStructSolidTet4NodeStab::projectStress(int varindex, double* outval)
{
  return;
}



void LagrangeElem3DStructSolidTet4NodeStab::projectStrain(int vartype, int varindex, double* outval)
{
  return;
}



void LagrangeElem3DStructSolidTet4NodeStab::projectIntVar(int index, double* outval)
{
  return;
}


int LagrangeElem3DStructSolidTet4NodeStab::calcOutput(double u1, double v1)
{
  return 0;
}



void LagrangeElem3DStructSolidTet4NodeStab::toPostprocess(int vartype, int varindex, int type, SparseMatrixXd&  coeffMat, VectorXd& rhsVec)
{
  return;
}


