
#include "LagrangeElem3DBbarFbar.h"

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



LagrangeElem3DBbarFbar::LagrangeElem3DBbarFbar()
{
  if (debug) cout << " constructor LagrangeElem3DBbarFbar\n\n";

  ndim   = 3;
  degree = 1;
  npElem = 8;
  nlbf   = 8;
  ndof   = 3;
  nsize  = nlbf*ndof;
}



LagrangeElem3DBbarFbar::~LagrangeElem3DBbarFbar()
{
  if (debug) cout << " destructor LagrangeElem3DBbarFbar\n\n";
}


void LagrangeElem3DBbarFbar::prepareElemData()
{
  LagrangeElement::prepareElemData();

  return;
}


void LagrangeElem3DBbarFbar::prepareElemData2()
{
  return;
}


int LagrangeElem3DBbarFbar::calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal, bool firstIter)
{
  if(finite)
    calcStiffnessAndResidualFS(Klocal, Flocal);
  else
    calcStiffnessAndResidualSS(Klocal, Flocal);

  return 0;
}




int LagrangeElem3DBbarFbar::calcStiffnessAndResidualSS(MatrixXd& Klocal, VectorXd& Flocal)
{
  double  F[9], Fc[9], detF, detFc, trFc, fact, fact1, fact2, fact3, dvol, dvol0;
  double  Jac, totvol=0.0;
  double  bb1, bb2, bb3, bb4, bb5, cc1, cc2, cc3, cc4, cc5;

  VectorXd  N(nlbf), dN_dx(nlbf), dN_dy(nlbf), dN_dz(nlbf);

  double  stre[6], cc[6][6], param[3], force[3], bforce[3], sig[3], acceCur[3];

  int   err,  isw,  count,  count1, index, ll = 0, ii, jj, kk, gp;
  int   TI, TIp1, TIp2, TJ, TJp1, TJp2;
  int   ind1, ind2;

  double rho0 = elmDat[5] ;
  double rho  = rho0 ;
  bforce[0]   = elmDat[6]*timeFunction[0].prop ;
  bforce[1]   = elmDat[7]*timeFunction[0].prop ;
  bforce[2]   = elmDat[8]*timeFunction[0].prop ;
  double  dt  = mpapTime.dt;
  double  af  = SolnData->td(2);
  double  acceFact1 = SolnData->td(5);
  double  acceFact2 = acceFact1;

  double xNode[8], yNode[8], zNode[8], xx, yy, zz;

  for(ii=0;ii<npElem;ii++)
  {
    xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
    yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    zNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][2];
  }

  //  compute determinant tr(F) in element centre 

  vector<double>  gausspoints1, gausspoints2, gausspoints3, gaussweights;

  getGaussPoints1D(1, gausspoints1, gaussweights);

  param[0] = gausspoints1[0];
  param[1] = gausspoints1[0];
  param[2] = gausspoints1[0];

  GeomData->computeBasisFunctions3D(0, 2, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);

  Fc[0] = computeValueCur(0, dN_dx) + 1.0;
  Fc[3] = computeValueCur(0, dN_dy);
  Fc[6] = computeValueCur(0, dN_dz);

  Fc[1] = computeValueCur(1, dN_dx);
  Fc[4] = computeValueCur(1, dN_dy) + 1.0;
  Fc[7] = computeValueCur(1, dN_dz);

  Fc[2] = computeValueCur(2, dN_dx);
  Fc[5] = computeValueCur(2, dN_dy);
  Fc[8] = computeValueCur(2, dN_dz) + 1.0;

  trFc  = Fc[0] + Fc[4] + Fc[8];
  detFc = Fc[0]*(Fc[4]*Fc[8]-Fc[5]*Fc[7]) - Fc[3]*(Fc[1]*Fc[8]-Fc[2]*Fc[7]) + Fc[6]*(Fc[1]*Fc[5]-Fc[2]*Fc[4]);

  VectorXd   Wmat1(nlbf), Wmat2(nlbf), Wmat3(nlbf);
  for(jj=0;jj<nlbf;jj++)
  {
    Wmat1(jj) = dN_dx[jj] ;
    Wmat2(jj) = dN_dy[jj] ;
    Wmat3(jj) = dN_dz[jj] ;
  }

  MatrixXd  Dmat(6,6), Bbar1(3,6), Bbar2(6,3), Ktemp(3,3), bc(3,6);

  getGaussPointsHex(nGP, gausspoints1, gausspoints2, gausspoints3, gaussweights);

  if(Klocal.rows() != nsize)
  {
    Klocal.resize(nsize, nsize);
    Flocal.resize(nsize);
  }
  Klocal.setZero();
  Flocal.setZero();

  //cout << " AAAAAAAAAAAAAAA " << endl;

  count = 1;   ll = 0;   err = 0;   isw = 3;
  totvol = 0.0;
  for(gp=0; gp<nGP; gp++)
  {
      param[2] = gausspoints3[gp];
      param[1] = gausspoints2[gp];
      param[0] = gausspoints1[gp];

        GeomData->computeBasisFunctions3D(0, 2, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);

        dvol0 = Jac * gaussweights[gp];
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

        acceCur[0] = computeValueDotDotCur(0, N);
        acceCur[1] = computeValueDotDotCur(1, N);
        acceCur[2] = computeValueDotDotCur(2, N);

        detF = F[0]*(F[4]*F[8]-F[5]*F[7]) - F[3]*(F[1]*F[8]-F[2]*F[7]) + F[6]*(F[1]*F[5]-F[2]*F[4]);

        if(detF < 0.0)   return 1;

        // ADJUST F33 fOR 2D PROBLEMS BASED ON THE ASSUMPTIONS OF PLANE STRESS/PLANE STRAIN/AXISYMMETRIC

        fact = (trFc - F[0] - F[4] - F[8])/3.0;

        F[0] += fact ;
        F[4] += fact ;
        F[8] += fact ;

        matlib3d_(matDat, F, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &isw, &err, &count, NULL);
        count++;
        ll += nivGP;

        if(err !=0)    return 1;

        for(ii=0;ii<6;ii++)
        {
          stre[ii] *= dvol;
          for(jj=0;jj<6;jj++)
          {
            cc[ii][jj] *= dvol;
            Dmat(ii,jj) = cc[ii][jj] ;
          }
        }

        //==============================================
        // CALCULATE TANGENT STIFFNESS
        //==============================================

        force[0] = bforce[0]*rho0;
        force[1] = bforce[1]*rho0;
        force[2] = bforce[2]*rho0;

        force[0] -= rho0*acceCur[0];
        force[1] -= rho0*acceCur[1];
        force[2] -= rho0*acceCur[2];

        for(ii=0;ii<nlbf;ii++)
        {
          bb1 = dN_dx[ii];
          bb2 = dN_dy[ii];
          bb3 = dN_dz[ii];
          bb4 = N[ii]*dvol0;

          fact1 = (Wmat1(ii) - bb1)/3.0;
          fact2 = (Wmat2(ii) - bb2)/3.0;
          fact3 = (Wmat3(ii) - bb3)/3.0;

          Bbar1(0,0) = bb1+fact1;
          Bbar1(0,1) = fact1;
          Bbar1(0,2) = fact1;
          Bbar1(0,3) = bb2;
          Bbar1(0,4) = 0.0;
          Bbar1(0,5) = bb3;

          Bbar1(1,0) = fact2;
          Bbar1(1,1) = bb2+fact2;
          Bbar1(1,2) = fact2;
          Bbar1(1,3) = bb1;
          Bbar1(1,4) = bb3;
          Bbar1(1,5) = 0.0;

          Bbar1(2,0) = fact3;
          Bbar1(2,1) = fact3;
          Bbar1(2,2) = bb3+fact3;
          Bbar1(2,3) = 0.0;
          Bbar1(2,4) = bb2;
          Bbar1(2,5) = bb1;

          bc = Bbar1*Dmat;

          TI   = 3*ii;
          TIp1 = TI+1;
          TIp2 = TI+2;

          Flocal(TI)   += bb4*force[0];
          Flocal(TIp1) += bb4*force[1];
          Flocal(TIp2) += bb4*force[2];

          Flocal(TI)    -= (Bbar1(0,0)*stre[0] + Bbar1(0,1)*stre[1] + Bbar1(0,2)*stre[2] + Bbar1(0,3)*stre[3] + Bbar1(0,4)*stre[4] + Bbar1(0,5)*stre[5]) ;
          Flocal(TIp1)  -= (Bbar1(1,0)*stre[0] + Bbar1(1,1)*stre[1] + Bbar1(1,2)*stre[2] + Bbar1(1,3)*stre[3] + Bbar1(1,4)*stre[4] + Bbar1(1,5)*stre[5]) ;
          Flocal(TIp2)  -= (Bbar1(2,0)*stre[0] + Bbar1(2,1)*stre[1] + Bbar1(2,2)*stre[2] + Bbar1(2,3)*stre[3] + Bbar1(2,4)*stre[4] + Bbar1(2,5)*stre[5]) ;

          for(jj=0; jj<nlbf; jj++)
          {
            cc1 = dN_dx[jj];
            cc2 = dN_dy[jj];
            cc3 = dN_dz[jj];
            cc4 = N[jj];

            fact1 = (Wmat1(jj) - cc1)/3.0;
            fact2 = (Wmat2(jj) - cc2)/3.0;
            fact3 = (Wmat3(jj) - cc3)/3.0;

            Bbar2(0,0) = cc1+fact1;
            Bbar2(1,0) = fact1;
            Bbar2(2,0) = fact1;
            Bbar2(3,0) = cc2;
            Bbar2(4,0) = 0.0;
            Bbar2(5,0) = cc3;

            Bbar2(0,1) = fact2;
            Bbar2(1,1) = cc2+fact2;
            Bbar2(2,1) = fact2;
            Bbar2(3,1) = cc1;
            Bbar2(4,1) = cc3;
            Bbar2(5,1) = 0.0;

            Bbar2(0,2) = fact3;
            Bbar2(1,2) = fact3;
            Bbar2(2,2) = cc3+fact3;
            Bbar2(3,2) = 0.0;
            Bbar2(4,2) = cc2;
            Bbar2(5,2) = cc1;

            TJ   = 3*jj;
            TJp1 = TJ+1;
            TJp2 = TJ+2;

            acceFact2 = acceFact1*cc4*rho0;

            fact  = bb4*acceFact2;

            Klocal(TI,   TJ)    +=  fact;
            Klocal(TIp1, TJp1)  +=  fact;
            Klocal(TIp2, TJp2)  +=  fact;

            Ktemp = (af*bc)*Bbar2;

            Klocal(TI,   TJ)    += Ktemp(0,0) ;
            Klocal(TI,   TJp1)  += Ktemp(0,1) ;
            Klocal(TI,   TJp2)  += Ktemp(0,2) ;

            Klocal(TIp1, TJ)    += Ktemp(1,0) ;
            Klocal(TIp1, TJp1)  += Ktemp(1,1) ;
            Klocal(TIp1, TJp2)  += Ktemp(1,2) ;

            Klocal(TIp2, TJ)    += Ktemp(2,0) ;
            Klocal(TIp2, TJp1)  += Ktemp(2,1) ;
            Klocal(TIp2, TJp2)  += Ktemp(2,2) ;
          }
        }
  }//gp

  return 0;
}




int LagrangeElem3DBbarFbar::calcStiffnessAndResidualFS(MatrixXd& Klocal, VectorXd& Flocal)
{
  double  F[9], Fc[9], detF, detFc, trFc, fact, fact1, fact2, fact3, dvol, dvol0;
  double  Jac, totvol=0.0, JacY, r1d3 = 1.0/3.0, r2d3=2.0/3.0;
  double  bb1, bb2, bb3, bb4, bb5, cc1, cc2, cc3, cc4, cc5;

  VectorXd  N(nlbf), dN_dx(nlbf), dN_dy(nlbf), dN_dz(nlbf);

  double  stre[6], cc[6][6], bc[3][6], param[3];
  double  force[3], bforce[3], sig[3], acceCur[3], cch[6];

  int   err,  isw,  count,  count1, index, ll, ii, jj, kk, gp;
  int   TI, TIp1, TIp2, TJ, TJp1, TJp2;
  int   ind1, ind2;

  double rho0 = elmDat[5] ;
  double rho  = rho0 ;
  bforce[0]   = elmDat[6]*timeFunction[0].prop ;
  bforce[1]   = elmDat[7]*timeFunction[0].prop ;
  bforce[2]   = elmDat[8]*timeFunction[0].prop ;
  double  dt  = mpapTime.dt;
  double  af  = SolnData->td(2);
  double  acceFact1 = SolnData->td(5);
  double  acceFact2 = acceFact1;


  double xNode[8], yNode[8], zNode[8], xx, yy, zz;

  for(ii=0;ii<npElem;ii++)
  {
    xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
    yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    zNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][2];
  }

  //totvol = (xNode[1]-xNode[0])*(yNode[2]-yNode[0])*(zNode[5]-zNode[0]);

  //  compute deformation gradient at the element center

  vector<double>  gausspoints1, gausspoints2, gausspoints3, gaussweights;

  getGaussPoints1D(1, gausspoints1, gaussweights);

  param[0] = gausspoints1[0];
  param[1] = gausspoints1[0];
  param[2] = gausspoints1[0];

  GeomData->computeBasisFunctions3D(0, 2, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);

  Fc[0] = computeValueCur(0, dN_dx) + 1.0;
  Fc[3] = computeValueCur(0, dN_dy);
  Fc[6] = computeValueCur(0, dN_dz);

  Fc[1] = computeValueCur(1, dN_dx);
  Fc[4] = computeValueCur(1, dN_dy) + 1.0;
  Fc[7] = computeValueCur(1, dN_dz);

  Fc[2] = computeValueCur(2, dN_dx);
  Fc[5] = computeValueCur(2, dN_dy);
  Fc[8] = computeValueCur(2, dN_dz) + 1.0;

  trFc  = Fc[0] + Fc[4] + Fc[8];
  detFc = Fc[0]*(Fc[4]*Fc[8] - Fc[5]*Fc[7]) - Fc[3]*(Fc[1]*Fc[8] - Fc[2]*Fc[7]) + Fc[6]*(Fc[1]*Fc[5] - Fc[2]*Fc[4]);

  if(finite)
    GeomData->computeBasisFunctions3D(1, 2, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);

  VectorXd   Wmat1(nlbf), Wmat2(nlbf), Wmat3(nlbf);
  for(jj=0;jj<nlbf;jj++)
  {
    Wmat1(jj) = dN_dx[jj] ;
    Wmat2(jj) = dN_dy[jj] ;
    Wmat3(jj) = dN_dz[jj] ;
  }

  //MatrixXd  Dmat(6,6), Bbar1(3,6), Bbar2(6,3), Ktemp(3,3);
  getGaussPointsHex(nGP, gausspoints1, gausspoints2, gausspoints3, gaussweights);

  if(Klocal.rows() != nsize)
  {
    Klocal.resize(nsize, nsize);
    Flocal.resize(nsize);
  }
  Klocal.setZero();
  Flocal.setZero();


  count = 1;   ll = 0;   err = 0;   isw = 3;
  totvol = 0.0;
  for(gp=0; gp<nGP; gp++)
  {
      param[2] = gausspoints3[gp];
      param[1] = gausspoints2[gp];
      param[0] = gausspoints1[gp];

        GeomData->computeBasisFunctions3D(0, 2, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);

        dvol0 = Jac*gaussweights[gp];
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

        acceCur[0] = computeValueDotDotCur(0, N);
        acceCur[1] = computeValueDotDotCur(1, N);
        acceCur[2] = computeValueDotDotCur(2, N);

        detF = F[0]*(F[4]*F[8] - F[5]*F[7]) - F[3]*(F[1]*F[8] - F[2]*F[7]) + F[6]*(F[1]*F[5] - F[2]*F[4]);

        if(detF < 0.0)   return 1;

        //  replace detF (tr(F)) with detF (tr(F)) from element centre

        if(finite)
        {
          fact = pow(detFc/detF, 1.0/3.0);

          for(ii=0; ii<9; ii++)
            F[ii] *= fact ;
        }
        else
        {
          fact = (trFc - F[0] - F[4] - F[8])/3.0;

          F[0] += fact ;
          F[4] += fact ;
          F[8] += fact ;
        }

        matlib3d_(matDat, F, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &isw, &err, &count, NULL);
        count++;
        ll += nivGP;

        if(err !=0)    return 1;

        if(finite)
        {
          GeomData->computeBasisFunctions3D(1, 2, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);

          dvol = Jac*gaussweights[gp];
        }

        for(ii=0;ii<6;ii++)
        {
          cch[ii]  =  r1d3 * (cc[ii][0] + cc[ii][1] + cc[ii][2]);
        }

        //==============================================
        // CALCULATE TANGENT STIFFNESS
        //==============================================

        force[0] = bforce[0]*rho0;
        force[1] = bforce[1]*rho0;
        force[2] = bforce[2]*rho0;

        force[0] -= rho0*acceCur[0];
        force[1] -= rho0*acceCur[1];
        force[2] -= rho0*acceCur[2];

        for(ii=0;ii<nlbf;ii++)
        {
          bb1 = dN_dx[ii]*dvol;
          bb2 = dN_dy[ii]*dvol;
          bb3 = dN_dz[ii]*dvol;
          bb4 = N[ii]*dvol;
          bb5 = N[ii]*dvol0;

          TI   = 3*ii;
          TIp1 = TI+1;
          TIp2 = TI+2;

          for(kk=0;kk<6;kk++)
          {
            bc[0][kk] = bb1*cc[0][kk] + bb2*cc[3][kk] + bb3*cc[5][kk];
            bc[1][kk] = bb1*cc[3][kk] + bb2*cc[1][kk] + bb3*cc[4][kk];
            bc[2][kk] = bb1*cc[5][kk] + bb2*cc[4][kk] + bb3*cc[2][kk];
          }

          sig[0] = bb1*stre[0] + bb2*stre[3] + bb3*stre[5] ;
          sig[1] = bb1*stre[3] + bb2*stre[1] + bb3*stre[4] ;
          sig[2] = bb1*stre[5] + bb2*stre[4] + bb3*stre[2] ;

          Flocal(TI)   += (bb5*force[0] - sig[0]) ;
          Flocal(TIp1) += (bb5*force[1] - sig[1]) ;
          Flocal(TIp2) += (bb5*force[2] - sig[2]) ;

          // derivative with respect to shp(centroid)

          fact1 = bb1 * cch[0] + bb2 * cch[3] + bb3 * cch[5] ;
          fact2 = bb1 * cch[3] + bb2 * cch[1] + bb3 * cch[4] ; 
          fact3 = bb1 * cch[5] + bb2 * cch[4] + bb3 * cch[2] ; 

          for(jj=0; jj<nlbf; jj++)
          {
            cc1 = dN_dx[jj];
            cc2 = dN_dy[jj];
            cc3 = dN_dz[jj];
            cc4 = N[jj];

            TJ   = 3*jj;
            TJp1 = TJ+1;
            TJp2 = TJ+2;

            acceFact2 = acceFact1*cc4*rho0;

            fact  = bb5*acceFact2;
            fact += af*(sig[0]*cc1 + sig[1]*cc2 + sig[2]*cc3);

            Klocal(TI,   TJ)    +=  fact;
            Klocal(TIp1, TJp1)  +=  fact;
            Klocal(TIp2, TJp2)  +=  fact;

            Klocal(TI,   TJ)    +=  af*(bc[0][0] * cc1 + bc[0][3] * cc2 + bc[0][5] * cc3) ;
            Klocal(TI,   TJp1)  +=  af*(bc[0][1] * cc2 + bc[0][3] * cc1 + bc[0][4] * cc3) ;
            Klocal(TI,   TJp2)  +=  af*(bc[0][2] * cc3 + bc[0][4] * cc2 + bc[0][5] * cc1) ;

            Klocal(TIp1, TJ)    +=  af*(bc[1][0] * cc1 + bc[1][3] * cc2 + bc[1][5] * cc3) ;
            Klocal(TIp1, TJp1)  +=  af*(bc[1][1] * cc2 + bc[1][3] * cc1 + bc[1][4] * cc3) ;
            Klocal(TIp1, TJp2)  +=  af*(bc[1][2] * cc3 + bc[1][4] * cc2 + bc[1][5] * cc1) ;

            Klocal(TIp2, TJ)    +=  af*(bc[2][0] * cc1 + bc[2][3] * cc2 + bc[2][5] * cc3) ;
            Klocal(TIp2, TJp1)  +=  af*(bc[2][1] * cc2 + bc[2][3] * cc1 + bc[2][4] * cc3) ;
            Klocal(TIp2, TJp2)  +=  af*(bc[2][2] * cc3 + bc[2][4] * cc2 + bc[2][5] * cc1) ;

            cc1 = af*(Wmat1(jj) - cc1);
            cc2 = af*(Wmat2(jj) - cc2);
            cc3 = af*(Wmat3(jj) - cc3);

            Klocal(TI,   TJ)    +=  fact1*cc1;
            Klocal(TI,   TJp1)  +=  fact1*cc2;
            Klocal(TI,   TJp2)  +=  fact1*cc3;

            Klocal(TIp1, TJ)    +=  fact2*cc1;
            Klocal(TIp1, TJp1)  +=  fact2*cc2;
            Klocal(TIp1, TJp2)  +=  fact2*cc3;

            Klocal(TIp2, TJ)    +=  fact3*cc1;
            Klocal(TIp2, TJp1)  +=  fact3*cc2;
            Klocal(TIp2, TJp2)  +=  fact3*cc3;

          }
        }
  }//gp

  return 0;
}



int LagrangeElem3DBbarFbar::calcLoadVector()
{
  return 0;
}



int LagrangeElem3DBbarFbar::calcInternalForces()
{
  return 0;
}



void LagrangeElem3DBbarFbar::discreteContourplot(int vartype, int varindex, int index, int nCol, double umin, double umax)
{
  return;
}


void LagrangeElem3DBbarFbar::projectToKnots(bool extrapolateFlag, int vartype, int varindex, int index)
{
  return;
}


void LagrangeElem3DBbarFbar::projectStress(int varindex, double* outval)
{
  return;
}



void LagrangeElem3DBbarFbar::projectStrain(int vartype, int varindex, double* outval)
{
  return;
}



void LagrangeElem3DBbarFbar::projectIntVar(int index, double* outval)
{
  return;
}


int LagrangeElem3DBbarFbar::calcOutput(double u1, double v1)
{
  return 0;
}



void LagrangeElem3DBbarFbar::toPostprocess(int vartype, int varindex, int type, SparseMatrixXd&  coeffMat, VectorXd& rhsVec)
{

  return;
}



void  LagrangeElem3DBbarFbar::computeEnergy(int index1, int index2, VectorXd& energy)
{
  cerr << " 'LagrangeElem3DBbarFbar::computeEnergy' is not complete. Need to correct it..."  << endl;

  // compute total energy for structural dynamic problems
  ///////////////////////////////////////////////////////////

  //elmDat = &(SolnData->ElemProp.data[0]);
  //matDat = &(SolidSolnData->MatlProp.data[0]);

  double  F[9], detF=0.0, F33, fact, dvol, dvol0, Jac, JacMult;
  double  posi[3], disp[3], velo[3];
  double  stress[6], strain[6], cc[6][6], param[3], bforce[3];

  VectorXd  N(nlbf), dN_dx(nlbf), dN_dy(nlbf);

  int   err,  isw,  count,  count1, index, ll = 0, ii, jj, kk, gp;
  int   TI, TIp1, TIp2, TJ, TJp1, TJp2;
  int   ind1, ind2;

  double  rho0  = elmDat[5] ;
  double  rho   = rho0 ;
  double  dt    = mpapTime.dt;

  vector<double>  gausspoints1, gausspoints2, gausspoints3, gaussweights;

  getGaussPointsHex(nGP, gausspoints1, gausspoints2, gausspoints3, gaussweights);

  energy.setZero();

  count = 1;   ll = 0;   err = 0;   isw = 3;
  for(gp=0; gp<nGP; gp++)
  {
      param[2] = gausspoints3[gp];
      param[1] = gausspoints2[gp];
      param[0] = gausspoints1[gp];

        GeomData->computeBasisFunctions2D(0, 2, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

        dvol0 = Jac*gaussweights[gp];
        dvol  = dvol0;

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
  }//gp

  printf("\n\n");

  return;
}





