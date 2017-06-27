
#include "LagrangeElem3DStructSolidHex8NodeStab.h"

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



LagrangeElem3DStructSolidHex8NodeStab::LagrangeElem3DStructSolidHex8NodeStab()
{
  if (debug) cout << " constructor LagrangeElem3DStructSolidHex8NodeStab\n\n";

  ndim   = 3;
  degree = 1;
  npElem = 8;
  nlbf   = 8;
  ndof   = 4;
  nsize  = nlbf*ndof;

}

LagrangeElem3DStructSolidHex8NodeStab::~LagrangeElem3DStructSolidHex8NodeStab()
{
  if (debug) cout << " destructor LagrangeElem3DStructSolidHex8NodeStab\n\n";
}


void LagrangeElem3DStructSolidHex8NodeStab::prepareElemData()
{
  LagrangeElement::prepareElemData();

  int ii, jj;

  ndim   = 3;
  degree = 1;
  npElem = 8;
  nlbf   = 8;
  ndof   = 4;
  nsize  = nlbf*ndof;

  // set the element property variables

  elmDat = &(SolnData->ElemProp[elmType].data[0]);
  matDat = &(SolnData->MatlProp[matType].data[0]);

  int nGP1  = (int) elmDat[0] ;
  int nGP2  = nGP1;
  nGP   = nGP1 * nGP2;
  finiteInt  = (int) elmDat[2] ;
  matId      = SolnData->MatlProp[matType].id + 1;
  finite     = (finiteInt >= 1) ;

  return;
}



void LagrangeElem3DStructSolidHex8NodeStab::prepareElemData2()
{
  cout << " LagrangeElem3DStructSolidHex8NodeStab::prepareElemData2() " << endl;
  return;
}


int LagrangeElem3DStructSolidHex8NodeStab::calcLoadVector()
{
  return 0;
}


int LagrangeElem3DStructSolidHex8NodeStab::calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal)
{
//   char fct[] = "LagrangeElem3DStructSolidHex8NodeStab::calcStiffnessAndResidual";
//   computerTime.go(fct);

  double  F[9], Fc[9], detF, detFc, trFc, fact, fact1, fact2, fact3, dvol, dvol0;
  double  Jac, JacTemp, totvol=0.0, JacY, r1d3 = 1.0/3.0, r2d3=2.0/3.0;
  double  bb1, bb2, bb3, bb4, bb5, cc1, cc2, cc3, cc4, cc5;
  double  pres, pbar, dp[3], tau, h, h2, volstr, rStab[3];

  VectorXd  N(nlbf), dN_dx(nlbf), dN_dy(nlbf), dN_dz(nlbf);
  MatrixXd  Dj(3,4);

  double  stre[6], strdev[6], cc[6][6], cctmp[6][6], bc[3][6], param[3], Idev[6][6];
  double  force[3], bforce[3], sig[3], acceCur[3];

  Idev3D(Idev);

  int   err,  isw,  count,  count1, index, ll, ii, jj, kk, mm, gp1, gp2, gp3;
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

  double  dt = mpapTime.dt;
  double  af = SolnData->td(2);
  double  acceFact1 = SolnData->td(5);
  double  acceFact2 = acceFact1;


  double xNode[8], yNode[8], zNode[8], xx, yy, zz;

  for(ii=0;ii<npElem;ii++)
  {
    xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
    yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    zNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][2];
  }

  totvol = (xNode[1]-xNode[0])*(yNode[2]-yNode[0])*(zNode[5]-zNode[0]);

  //h = pow(6.0*totvol/PI, 1.0/3.0);
  
  h = xNode[1]-xNode[0];
  //h = sqrt(3.0*h*h);
  
  double  alpha = 1.0/100.0;

  //tau = elmDat[8]*alpha*h*h/sqrt(mu);
  tau = elmDat[8]*alpha*h*h/(mu*dt);

  //tau = elmDat[8]*dt*0.3/0.9;

  vector<double>  gausspoints, gaussweights;

  //MatrixXd  Dmat(6,6), Bbar1(3,6), Bbar2(6,3), Ktemp(3,3);

  count = 1;   ll = 0;   err = 0;   isw = 3;

  nGP = (int) elmDat[0] ;

  getGaussPoints1D(nGP, gausspoints, gaussweights);

    if(Klocal.rows() != nsize)
    {
      Klocal.resize(nsize, nsize);
      Flocal.resize(nsize);
    }
    Klocal.setZero();
    Flocal.setZero();


  totvol = 0.0;
  for(gp3=0;gp3<nGP;gp3++)
  {
    param[2] = gausspoints[gp3];
  for(gp2=0;gp2<nGP;gp2++)
  {
    param[1] = gausspoints[gp2];

    JacY = gaussweights[gp3] * gaussweights[gp2];

  for(gp1=0;gp1<nGP;gp1++)
  {
        param[0] = gausspoints[gp1];

        GeomData->computeBasisFunctions3D(0, 2, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);

        JacTemp = JacY * gaussweights[gp1];

        dvol0 = Jac * JacTemp;
        dvol = dvol0;
        totvol += dvol;

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

        pres  = computeValueCur(3, N);

        acceCur[0] = computeValueDotDotCur(0, N);
        acceCur[1] = computeValueDotDotCur(1, N);
        acceCur[2] = computeValueDotDotCur(2, N);

        detF = F[0]*(F[4]*F[8] - F[5]*F[7]) - F[3]*(F[1]*F[8] - F[2]*F[7]) + F[6]*(F[1]*F[5] - F[2]*F[4]);

        if(detF < 0.0)   return 1;

        //  replace detF (tr(F)) with detF (tr(F)) from element centre

        matlib3d_(matDat, F, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &isw, &err, &count, NULL);
        count++;
        ll += nivGP;

        if(err !=0)    return 1;

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
          GeomData->computeBasisFunctions3D(1, 2, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);

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
          Flocal(TIp3) -= (bb5*volstr);

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

              Klocal(TIp3, TJ)    +=  af*(bb4 * cc1) ;
              Klocal(TIp3, TJp1)  +=  af*(bb4 * cc2) ;
              Klocal(TIp3, TJp2)  +=  af*(bb4 * cc3);
              Klocal(TIp3, TJp3)  -=  af*(bb5 * cc4)*eps;

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
  }//gp1
  }//gp2
  }//gp3

  return 0;
}






int LagrangeElem3DStructSolidHex8NodeStab::calcInternalForces()
{
  return 0;
}



void LagrangeElem3DStructSolidHex8NodeStab::discreteContourplot(int vartype, int varindex, int index, int nCol, double umin, double umax)
{
  return;
}


void LagrangeElem3DStructSolidHex8NodeStab::projectToKnots(bool extrapolateFlag, int vartype, int varindex, int index)
{
  return;
}


void LagrangeElem3DStructSolidHex8NodeStab::projectStress(int varindex, double* outval)
{
  return;
}



void LagrangeElem3DStructSolidHex8NodeStab::projectStrain(int vartype, int varindex, double* outval)
{
  return;
}



void LagrangeElem3DStructSolidHex8NodeStab::projectIntVar(int index, double* outval)
{
  return;
}


int LagrangeElem3DStructSolidHex8NodeStab::calcOutput(double u1, double v1)
{
  return 0;
}



void LagrangeElem3DStructSolidHex8NodeStab::toPostprocess(int vartype, int varindex, int type, SparseMatrixXd&  coeffMat, VectorXd& rhsVec)
{
  return;
}


void  LagrangeElem3DStructSolidHex8NodeStab::computeEnergy(int index1, int index2, VectorXd& energy)
{
  // compute total energy for structural dynamic problems
  ///////////////////////////////////////////////////////////

  //elmDat = &(SolnData->ElemProp.data[0]);
  //matDat = &(SolidSolnData->MatlProp.data[0]);

  double  F[4], detF=0.0, F33, fact, dvol, dvol0, Jac, dt, JacMult;
  double  posi[2], disp[2], velo[2];
  double  stress[4], strain[4], cc[4][4], param[2], bforce[2];

  VectorXd  N(nlbf), dN_dx(nlbf), dN_dy(nlbf);

  int   err,  isw,  count,  count1, index, ll = 0, ii, jj, gp1, gp2, TI, TIp1, TJ, TJp1;
  int   ind1, ind2, kk;
  
  double  rho  = elmDat[5] ;

  count = 1;   ll = 0;   err = 0;   isw = 3;
  dt = mpapTime.dt;

  //cout << " finite = " << finite << endl;

  int nGP = (int) elmDat[0] ;

  vector<double>  gausspoints1, gaussweights1;

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

        matlib2d_(matDat, F, &F33, stress, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);
        count++;
        ll += nivGP;

        strain[0] = F[0] - 1.0;
        strain[1] = F[3] - 1.0;
        strain[2] = F33  - 1.0;
        strain[3] = 0.5 * (F[1] + F[2]);

        printf(" %14.12f \t %14.12f \t %14.12f \t %14.12f \n ", strain[0], strain[1], strain[2], strain[3]);
        printf(" %14.12f \t %14.12f \t %14.12f \t %14.12f \n ", stress[0], stress[1], stress[2], stress[3]);
        //printf("  dvol =  %14.12f \n\n ", dvol);

        velo[0] = computeValueDot(0, N);
        velo[1] = computeValueDot(1, N);

        // kinetic energy
        energy[0] += (0.5*dvol0*rho*(velo[0]*velo[0]+velo[1]*velo[1]));
        energy[1] += (0.5*dvol*(strain[0]*stress[0]+strain[1]*stress[1]+strain[2]*stress[2]+strain[3]*stress[3]));

  }//gp1
  }//gp2

  printf("\n\n");

  return;
}




