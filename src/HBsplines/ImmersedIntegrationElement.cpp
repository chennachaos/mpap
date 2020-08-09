
#include "ImmersedIntegrationElement.h"
#include "BasisFunctionsLagrange.h"
#include "GeomDataLagrange.h"
#include "SolutionData.h"
#include "QuadratureUtil.h"
#include "myConstants.h"


ImmersedIntegrationElement::ImmersedIntegrationElement()
{
  id = itemcount++;
  DIM = 2;

  IS_ACTIVE = true;
}


ImmersedIntegrationElement::~ImmersedIntegrationElement()
{
  itemcount--;
}

void ImmersedIntegrationElement::prepareElemData()
{
  return;
}


void ImmersedIntegrationElement::initialiseDOFvalues()
{
  int  ii, jj, ind1;

  ind1 = pointNums.size();
  elemNums.resize(ind1);

  //cout << " ind1 " << ind1 << '\t' << DIM << endl;

  //printVector(posIndices);

  for(ii=0;ii<ind1;ii++)
  {
    for(jj=0;jj<DIM;jj++)
      posIndices.push_back(pointNums[ii]*DIM+jj);
  }

  //printVector(posIndices);

  nGP = ind1;

  if(nGP == 1)
  {
    gausspoints.push_back(1.0);
    gaussweights.push_back(1.0);
  }
  else
    getGaussPoints1D(nGP, gausspoints, gaussweights);

  return;
}


void ImmersedIntegrationElement::resetMatrixAndVector()
{
    Khorz.setZero();
    Kvert.setZero();
    Flocal.setZero();
    Flocal2.setZero();

   return;
}


void ImmersedIntegrationElement::reset()
{
  return;
}


void ImmersedIntegrationElement::assembleElementVector(int ind, bool flag, double* rhs)
{
  return;
}


void ImmersedIntegrationElement::assembleElementMatrix(int index, SparseMatrixXd& mtx)
{
  return;
}


void ImmersedIntegrationElement::assembleMatrixAndVector(int start, SparseMatrixXd& mtx, double* rhs)
{
  return;
}


void ImmersedIntegrationElement::computePointAtGP(int ind, myPoint& pt)
{
  int  ii, jj, nlb;
  nlb = pointNums.size();
  vector<double>  N(nlb);

  Lagrange_BasisFuns1D(nlb-1, gausspoints[ind], &N[0]);

  pt.setZero();
  for(ii=0;ii<nlb;ii++)
  {
    for(jj=0;jj<DIM;jj++)
      pt[jj] += (GeomDataLag->NodePosCur[pointNums[ii]][jj] * N[ii] );
  }

  return;
}


void ImmersedIntegrationElement::computeVelocity(const VectorXd& NN, myPoint&  velSpec)
{
  assert(NN.rows() == pointNums.size());

  int  ii, jj;
  velSpec.setZero();
  for(ii=0;ii<pointNums.size();ii++)
  {
    for(jj=0; jj<DIM; jj++)
      velSpec[jj] += (GeomDataLag->specValNew[pointNums[ii]][jj] * NN[ii] );
  }

  return;
}




void ImmersedIntegrationElement::computeVelocityCur(const VectorXd& NN, myPoint&  velSpec)
{
  assert(NN.rows() == pointNums.size());

  int  ii, jj;
  velSpec.setZero();
  for(ii=0;ii<pointNums.size();ii++)
  {
    for(jj=0; jj<DIM; jj++)
      velSpec[jj] += (GeomDataLag->specValCur[pointNums[ii]][jj] * NN[ii] );
  }

  return;
}




void ImmersedIntegrationElement::computeAcceleration(const VectorXd& NN, myPoint&  velSpec)
{
  assert(NN.rows() == pointNums.size());

  int  ii, jj;
  velSpec.setZero();
  for(ii=0;ii<pointNums.size();ii++)
  {
    for(jj=0; jj<DIM; jj++)
      velSpec[jj] += (GeomDataLag->acceNew[pointNums[ii]][jj] * NN[ii] );
  }

  return;
}




void ImmersedIntegrationElement::computeAccelerationCur(const VectorXd& NN, myPoint&  velSpec)
{
  assert(NN.rows() == pointNums.size());

  int  ii, jj;
  velSpec.setZero();
  for(ii=0;ii<pointNums.size();ii++)
  {
    for(jj=0; jj<DIM; jj++)
      velSpec[jj] += (GeomDataLag->acceCur[pointNums[ii]][jj] * NN[ii] );
  }

  return;
}




void ImmersedIntegrationElement::integrateForceAndMoment(VectorXd& vectemp, myPoint& centLoc)
{
  // to compute total force on the immersed rigid solid
  
  int  gp, ii, jj, kk, nlb = pointNums.size();
  double  dvol, detJ;
  vector<double>  N(nlb), dN(nlb), dN_dx(nlb), xx(nlb), yy(nlb), zz(nlb);
  double  lamLoc[3]={0.0, 0.0, 0.0}, rad, fact;
  double xc = centLoc[0];
  double yc = centLoc[1];
  double zc = centLoc[2];

  for(ii=0;ii<nlb;ii++)
  {
    xx[ii] = GeomDataLag->NodePosCur[pointNums[ii]][0];
    yy[ii] = GeomDataLag->NodePosCur[pointNums[ii]][1];
    zz[ii] = 0.0;
  }

  bool axsy = (GeomDataHBS->FluidProps[2] == 1);

  if(DIM == 3)
  {
    for(ii=0;ii<nlb;ii++)
      zz[ii] = GeomDataLag->NodePosCur[pointNums[ii]][2];

    axsy = false;
  }


  //cout << " cccccccccc " << endl;
  for(gp=0;gp<gausspoints.size();gp++)
  {
    computeLagrangeBFsLine3D(nlb-1, gausspoints[gp], &xx[0], &yy[0], &zz[0], &N[0], &dN[0], detJ);

    dvol  = gaussweights[gp] * detJ;

    if(axsy)
    {
      rad = 0.0;
      for(ii=0;ii<nlb;ii++)
        rad += N[ii] * xx[ii];

      dvol *= 2.0*PI*rad;
    }

    for(ii=0;ii<nlb;ii++)
    {
      //cout << ii << '\t' << N[ii] << '\t' << dvol << endl;
      kk = pointNums[ii]*DIM;

      for(jj=0;jj<DIM;jj++)
        lamLoc[jj] = SolnData->var3(kk+jj) ;

      //cout << lamLoc[0] << '\t' << lamLoc[1] << '\t' << lamLoc[2] << endl;

      fact = N[ii]*dvol;

      vectemp[0] += ( fact*lamLoc[0] );
      vectemp[1] += ( fact*lamLoc[1] );
      vectemp[2] += ( fact*lamLoc[2] );

      vectemp[3] += ( (yc*lamLoc[2]-zc*lamLoc[1]) * fact);
      vectemp[4] += ( (zc*lamLoc[0]-xc*lamLoc[2]) * fact);
      vectemp[5] += ( (xc*lamLoc[1]-yc*lamLoc[0]) * fact);
    }
  }

  return;
}




void ImmersedIntegrationElement::computeKhorzKvertRigid(MatrixXd& Kh, MatrixXd& Kv)
{
  int  gp, ii, nlb, size, kk;
  double  dvol, detJ, dvol1, dvol2, fact;
  nlb = pointNums.size();
  vector<double>  N(nlb), dN(nlb), dN_dx(nlb), xx(nlb), yy(nlb);
  double  dx, dy, xc, yc;

  for(ii=0;ii<nlb;ii++)
  {
    xx[ii] = GeomDataLag->NodePosCur[pointNums[ii]][0];
    yy[ii] = GeomDataLag->NodePosCur[pointNums[ii]][1];
  }

  bool axsy = (GeomDataHBS->FluidProps[2] == 1);
  
  size = 3; // in 2D

  Kh.resize(size, nlb*DIM);
  Kh.setZero();

  Kv.resize(nlb*DIM, size);
  Kv.setZero();

  for(gp=0;gp<gausspoints.size();gp++)
  {
    computeLagrangeBFsLine2D(nlb-1, gausspoints[gp], &xx[0], &yy[0], &N[0], &dN[0], detJ);

    //cout << " detJ " << detJ << '\t' << gaussweights[gp] << endl;

    dvol  = gaussweights[gp] * detJ;

    dvol1 = dvol;

    xc = yc = 0.0;
    for(ii=0;ii<nlb;ii++)
    {
      xc += N[ii] * xx[ii];
      yc += N[ii] * yy[ii];
    }

    if(axsy)
    {
      dvol1 = dvol * 2.0*PI*xc;

      //if(rad < 0.0)
      //{
        //dvol1 = dvol * -1.0;
        //dvol2 = dvol * -1.0;
      //}
      //
    }

    dvol2 = dvol;

    for(ii=0;ii<nlb;ii++)
    {
      fact = N[ii]*dvol1;

      kk = ii*DIM;
      Kh(0, kk)   += fact;
      Kh(1, kk+1) += fact;
      Kh(2, kk)   += ( -yc*fact);
      Kh(2, kk+1) += (  xc*fact);

      fact = N[ii]*dvol2;

      Kv(kk,   0) += fact;
      Kv(kk+1, 1) += fact;
      Kv(kk,   2) += ( -yc*fact);
      Kv(kk+1, 2) += (  xc*fact);
    }
  }

  return;
}





void ImmersedIntegrationElement::integrateForceFlexible(int ind1, int ind2, VectorXd& vectemp)
{
  // to compute force vector on the immersed flexible solid

  int  gp, ii, kk;
  int  nlbL = pointNums.size();
  int  nlbS = nlbL;

  vector<double>  N(nlbL), dN(nlbL), dN_dx(nlbL), xNode(nlbL), yNode(nlbL);
  double  xc, dvol, detJ, lamX, lamY, fact;

  for(ii=0;ii<nlbL;ii++)
  {
    xNode[ii] = GeomDataLag->NodePosCur[pointNums[ii]][0];
    yNode[ii] = GeomDataLag->NodePosCur[pointNums[ii]][1];
  }

  bool axsy = (GeomDataHBS->FluidProps[2] == 1);

  for(gp=0;gp<gausspoints.size();gp++)
  {
    computeLagrangeBFsLine2D(nlbL-1, gausspoints[gp], &xNode[0], &yNode[0], &N[0], &dN[0], detJ);

    dvol  = gaussweights[gp] * detJ;

    if(axsy)
    {
      xc = 0.0;
      for(ii=0; ii<nlbL; ii++)
        xc += N[ii] * xNode[ii];

      dvol *= 2.0*PI*xc;

      //if(xc < 0.0)
        //dvol *= -1.0;
    }

    lamX = 0.0; // force in X-direction
    lamY = 0.0; // force in Y-direction
    for(ii=0; ii<nlbL; ii++)
    {
      //cout << ii << '\t' << N[ii] << endl;
      kk = pointNums[ii]*DIM;

      lamX  += SolnData->var3(kk)   * N[ii] ;
      lamY  += SolnData->var3(kk+1) * N[ii] ;
    }

    for(ii=0; ii<nlbS; ii++)
    {
      fact = N[ii]*dvol;

      kk = pointNums[ii]*3;

      vectemp[kk]   += ( fact * lamX); // force in X-direction
      vectemp[kk+1] += ( fact * lamY); // force in Y-direction
    }
  }

  return;
}






void ImmersedIntegrationElement::computeKhorzKvertFlexible(int ind1, int ind2, MatrixXd& Kh, MatrixXd& Kv)
{
  int  gp, ii, jj, TI, TIp1, TIp2, TJ, TJp1;

  int  nlbL = pointNums.size();
  int  nlbS = nlbL;

  vector<double>  NL(nlbL), dN(nlbL), dN_dx(nlbL), xNode(nlbL), yNode(nlbL), NS(nlbS);
  double  dvol, detJ, dvol1, dvol2, fact1, fact2, xc;

  for(ii=0;ii<nlbL;ii++)
  {
    xNode[ii] = GeomDataLag->NodePosCur[pointNums[ii]][0];
    yNode[ii] = GeomDataLag->NodePosCur[pointNums[ii]][1];
  }

  bool axsy = (GeomDataHBS->FluidProps[2] == 1);
  
  Kh.resize(3*nlbS, DIM*nlbL);
  Kh.setZero();

  Kv.resize(DIM*nlbS, 3*nlbS);
  Kv.setZero();

  for(gp=0;gp<gausspoints.size();gp++)
  {
    computeLagrangeBFsLine2D(nlbL-1, gausspoints[gp], &xNode[0], &yNode[0], &NL[0], &dN[0], detJ);

    //cout << " detJ " << detJ << '\t' << gaussweights[gp] << endl;

    dvol  = gaussweights[gp] * detJ;
    
    dvol1 = dvol;

    xc = 0.0;
    for(ii=0;ii<nlbL;ii++)
    {
      xc += NL[ii] * xNode[ii];
    }

    if(axsy)
    {
      dvol1 = dvol * 2.0*PI*xc;

      //if(rad < 0.0)
      //{
        //dvol1 = dvol * -1.0;
        //dvol2 = dvol * -1.0;
      //}
      //
    }

    dvol2 = dvol1;

    if(nlbS == nlbL)
    {
      for(ii=0; ii<nlbL; ii++)
        NS[ii] = NL[ii];
    }
    else
    {
      NS[0] = 1.0;
    }

    for(ii=0; ii<nlbS; ii++)
    {
      fact1 = NS[ii]*dvol1;
      fact2 = NS[ii]*dvol2;

      TI   = ii*3;
      TIp1 = TI+1;
      TIp2 = TI+2;

      for(jj=0; jj<nlbL; jj++)
      {
        TJ   = jj*DIM;
        TJp1 = TJ+1;
        
        fact1 *= NL[jj];

        Kh(TI,   TJ)   += fact1;
        Kh(TIp1, TJp1) += fact1;
        //Kh(TIp2, TJ)   += ( -yc*fact1);
        //Kh(TIp2, TJp1) += (  xc*fact1);

        fact2 *= NL[jj];

        Kv(TJ,   TI)   += fact2;
        Kv(TJp1, TIp1) += fact2;
        //Kv(TJ,   TIp2) += ( -yc*fact2);
        //Kv(TJp1, TIp2) += (  xc*fact2);
      }
    }

  }

  return;
}








