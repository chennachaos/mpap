
#include "TreeNode.h"
#include "MpapTime.h"
#include "TimeFunction.h"
#include "Functions.h"
#include "SolutionData.h"
#include "BasisFunctionsBSpline.h"
#include "myDataIntegrateCutFEM.h"
#include "myPoly.h"


extern  MpapTime  mpapTime;
extern List<TimeFunction> timeFunction;



template<>
void TreeNode<2>::computeVelocity(myPoint& param, myPoint& vel)
{
    int ii, TI;

    double  b1, b2;

    VectorXd  NN(totnlbf), N;

    GeomData->computeBasisFunctions2D(knotBegin, knotIncr, param, NN);

    if(parent == NULL)
    {
      N = NN;
    }
    else
    {
      N = SubDivMat*NN;
    }

    vel.setZero();
    for(ii=0;ii<totnlbf2;ii++)
    {
      TI = ndof*GlobalBasisFuncs[ii];

      b1 = SolnData->var1(TI);
      b2 = SolnData->var1(TI+1);

      vel(0) += ( b1 * N(ii) );
      vel(1) += ( b2 * N(ii) );
    }

  return;
}



template<>
void TreeNode<2>::computeVelocityCur(myPoint& param, myPoint& vel)
{
    int ii, TI;

    double  b1, b2;

    VectorXd  NN(totnlbf), N;

    GeomData->computeBasisFunctions2D(knotBegin, knotIncr, param, NN);

    if(parent == NULL)
    {
      N = NN;
    }
    else
    {
      N = SubDivMat*NN;
    }

    vel.setZero();
    for(ii=0;ii<totnlbf2;ii++)
    {
      TI = ndof*GlobalBasisFuncs[ii];

      b1 = SolnData->var1Cur(TI);
      b2 = SolnData->var1Cur(TI+1);

      vel(0) += ( b1 * N(ii) );
      vel(1) += ( b2 * N(ii) );
    }

  return;
}



template<>
void TreeNode<3>::computeVelocity(myPoint& param, myPoint& vel)
{
    int ii, TI;

    double  b1, b2, b3;

    VectorXd  NN(totnlbf), N;

    GeomData->computeBasisFunctions3D(knotBegin, knotIncr, param, NN);

    if(parent == NULL)
    {
      N = NN;
    }
    else
    {
      N = SubDivMat*NN;
    }

    vel.setZero();
    for(ii=0;ii<totnlbf2;ii++)
    {
      TI = ndof*GlobalBasisFuncs[ii];

      b1 = SolnData->var1(TI);
      b2 = SolnData->var1(TI+1);
      b3 = SolnData->var1(TI+2);

      vel(0) += ( b1 * N(ii) );
      vel(1) += ( b2 * N(ii) );
      vel(2) += ( b3 * N(ii) );
    }

  return;
}



template<>
void TreeNode<3>::computeVelocityCur(myPoint& param, myPoint& vel)
{
    int ii, TI;

    double  b1, b2, b3;

    VectorXd  NN(totnlbf), N;

    GeomData->computeBasisFunctions3D(knotBegin, knotIncr, param, NN);

    if(parent == NULL)
    {
      N = NN;
    }
    else
    {
      N = SubDivMat*NN;
    }

    vel.setZero();
    for(ii=0;ii<totnlbf2;ii++)
    {
      TI = ndof*GlobalBasisFuncs[ii];

      b1 = SolnData->var1Cur(TI);
      b2 = SolnData->var1Cur(TI+1);
      b3 = SolnData->var1Cur(TI+2);

      vel(0) += ( b1 * N(ii) );
      vel(1) += ( b2 * N(ii) );
      vel(2) += ( b3 * N(ii) );
    }

  return;
}






template<>
void TreeNode<2>::computeVelocityGradient(myPoint& param, MatrixXd& velGrad)
{
    int ii, TI;

    double  b1, b2;

    VectorXd  NN(totnlbf), dNN_dx(totnlbf), dNN_dy(totnlbf);
    VectorXd  N, dN_dx, dN_dy;

    GeomData->computeBasisFunctions2D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy);

    if(parent == NULL)
    {
      N = NN;

      dN_dx = dNN_dx;
      dN_dy = dNN_dy;
    }
    else
    {
      N = SubDivMat*NN;

      dN_dx = SubDivMat*dNN_dx;
      dN_dy = SubDivMat*dNN_dy;
    }

    velGrad.setZero();
    for(ii=0;ii<totnlbf2;ii++)
    {
      TI = ndof*GlobalBasisFuncs[ii];

      b1 = SolnData->var1(TI);
      b2 = SolnData->var1(TI+1);

      velGrad(0,0) += ( b1 * dN_dx(ii) );
      velGrad(0,1) += ( b1 * dN_dy(ii) );

      velGrad(1,0) += ( b2 * dN_dx(ii) );
      velGrad(1,1) += ( b2 * dN_dy(ii) );
    }

  return;
}



template<>
void TreeNode<2>::computeVelocityGradientCur(myPoint& param, MatrixXd& velGrad)
{
    int ii, TI;

    double  b1, b2;

    VectorXd  NN(totnlbf), dNN_dx(totnlbf), dNN_dy(totnlbf);
    VectorXd  N, dN_dx, dN_dy;

    GeomData->computeBasisFunctions2D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy);

    if(parent == NULL)
    {
      N = NN;

      dN_dx = dNN_dx;
      dN_dy = dNN_dy;
    }
    else
    {
      N = SubDivMat*NN;

      dN_dx = SubDivMat*dNN_dx;
      dN_dy = SubDivMat*dNN_dy;
    }

    velGrad.setZero();
    for(ii=0;ii<totnlbf2;ii++)
    {
      TI = ndof*GlobalBasisFuncs[ii];

      b1 = SolnData->var1Cur(TI);
      b2 = SolnData->var1Cur(TI+1);

      velGrad(0,0) += ( b1 * dN_dx(ii) );
      velGrad(0,1) += ( b1 * dN_dy(ii) );

      velGrad(1,0) += ( b2 * dN_dx(ii) );
      velGrad(1,1) += ( b2 * dN_dy(ii) );
    }

  return;
}



template<>
void TreeNode<3>::computeVelocityGradient(myPoint& param, MatrixXd& velGrad)
{
    int ii, TI;

    double  b1, b2, b3;

    VectorXd  NN(totnlbf), dNN_dx(totnlbf), dNN_dy(totnlbf), dNN_dz(totnlbf);
    VectorXd  N, dN_dx, dN_dy, dN_dz;

    GeomData->computeBasisFunctions3D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy, dNN_dz);

    if(parent == NULL)
    {
      N = NN;

      dN_dx = dNN_dx;
      dN_dy = dNN_dy;
      dN_dz = dNN_dz;
    }
    else
    {
      N = SubDivMat*NN;

      dN_dx = SubDivMat*dNN_dx;
      dN_dy = SubDivMat*dNN_dy;
      dN_dz = SubDivMat*dNN_dz;
    }

    velGrad.setZero();
    for(ii=0;ii<totnlbf2;ii++)
    {
      TI = ndof*GlobalBasisFuncs[ii];

      b1 = SolnData->var1(TI);
      b2 = SolnData->var1(TI+1);
      b3 = SolnData->var1(TI+2);

      velGrad(0,0) += ( b1 * dN_dx(ii) );
      velGrad(0,1) += ( b1 * dN_dy(ii) );
      velGrad(0,2) += ( b1 * dN_dz(ii) );

      velGrad(1,0) += ( b2 * dN_dx(ii) );
      velGrad(1,1) += ( b2 * dN_dy(ii) );
      velGrad(1,2) += ( b2 * dN_dz(ii) );

      velGrad(2,0) += ( b3 * dN_dx(ii) );
      velGrad(2,1) += ( b3 * dN_dy(ii) );
      velGrad(2,2) += ( b3 * dN_dz(ii) );
    }

  return;
}


template<>
void TreeNode<3>::computeVelocityGradientCur(myPoint& param, MatrixXd& velGrad)
{
    int ii, TI;

    double  b1, b2, b3;

    VectorXd  NN(totnlbf), dNN_dx(totnlbf), dNN_dy(totnlbf), dNN_dz(totnlbf);
    VectorXd  N, dN_dx, dN_dy, dN_dz;

    GeomData->computeBasisFunctions3D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy, dNN_dz);

    if(parent == NULL)
    {
      N = NN;

      dN_dx = dNN_dx;
      dN_dy = dNN_dy;
      dN_dz = dNN_dz;
    }
    else
    {
      N = SubDivMat*NN;

      dN_dx = SubDivMat*dNN_dx;
      dN_dy = SubDivMat*dNN_dy;
      dN_dz = SubDivMat*dNN_dz;
    }

    velGrad.setZero();
    for(ii=0;ii<totnlbf2;ii++)
    {
      TI = ndof*GlobalBasisFuncs[ii];

      b1 = SolnData->var1Cur(TI);
      b2 = SolnData->var1Cur(TI+1);
      b3 = SolnData->var1Cur(TI+2);

      velGrad(0,0) += ( b1 * dN_dx(ii) );
      velGrad(0,1) += ( b1 * dN_dy(ii) );
      velGrad(0,2) += ( b1 * dN_dz(ii) );

      velGrad(1,0) += ( b2 * dN_dx(ii) );
      velGrad(1,1) += ( b2 * dN_dy(ii) );
      velGrad(1,2) += ( b2 * dN_dz(ii) );

      velGrad(2,0) += ( b3 * dN_dx(ii) );
      velGrad(2,1) += ( b3 * dN_dy(ii) );
      velGrad(2,2) += ( b3 * dN_dz(ii) );
    }

  return;
}



template<>
void TreeNode<1>::computeVelocityAndStress(myPoint& param, myPoint& vel, MatrixXd& stress, myPoint& acce)
{
  return;
}





template<>
void TreeNode<2>::computeVelocityAndStress(myPoint& param, myPoint& vel, MatrixXd& stress, myPoint& acce)
{
    int ii, TI;

    double  b1, b2, b3, b4, pres;

    VectorXd  NN(totnlbf), dNN_dx(totnlbf), dNN_dy(totnlbf), velDot(2);
    VectorXd  N, dN_dx, dN_dy;

    double  rho = elmDat[3];
    double  mu  = elmDat[4];

    GeomData->computeBasisFunctions2D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy);

    if(parent == NULL)
    {
      N = NN;

      dN_dx = dNN_dx;
      dN_dy = dNN_dy;
    }
    else
    {
      N = SubDivMat*NN;

      dN_dx = SubDivMat*dNN_dx;
      dN_dy = SubDivMat*dNN_dy;
    }

    vel.setZero();
    acce.setZero();
    stress.setZero();
    pres = 0.0;

    for(ii=0;ii<totnlbf2;ii++)
    {
      TI = ndof*GlobalBasisFuncs[ii];

      b1 = SolnData->var1(TI);
      b2 = SolnData->var1(TI+1);
      b4 = SolnData->var1(TI+2);

      vel(0) += ( b1 * N(ii) );
      vel(1) += ( b2 * N(ii) );

      stress(0,0) += ( b1 * dN_dx(ii) );
      stress(0,1) += ( b1 * dN_dy(ii) );

      stress(1,0) += ( b2 * dN_dx(ii) );
      stress(1,1) += ( b2 * dN_dy(ii) );

      pres   += ( b4 * N(ii) );

      b1  = SolnData->var1DotCur(TI);
      b2  = SolnData->var1DotCur(TI+1);

      acce(0) += ( b1 * N(ii) );
      acce(1) += ( b2 * N(ii) );
    }

    // this is pseudo-stress
    stress *= mu;
    //stress = mu*(stress+stress.transpose());

    stress(0,0) -= pres;
    stress(1,1) -= pres;

  return;
}



template<>
void TreeNode<2>::computeVelocityAndStressCur(myPoint& param, myPoint& vel, MatrixXd& stress, myPoint& acce)
{
    int ii, TI;

    double  b1, b2, b4, pres;

    VectorXd  NN(totnlbf), dNN_dx(totnlbf), dNN_dy(totnlbf);
    VectorXd  N, dN_dx, dN_dy;

    double  mu  = elmDat[4];

    GeomData->computeBasisFunctions2D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy);

    if(parent == NULL)
    {
      N = NN;

      dN_dx = dNN_dx;
      dN_dy = dNN_dy;
    }
    else
    {
      N = SubDivMat*NN;

      dN_dx = SubDivMat*dNN_dx;
      dN_dy = SubDivMat*dNN_dy;
    }

    vel.setZero();
    stress.setZero();
    pres = 0.0;

    for(ii=0;ii<totnlbf2;ii++)
    {
      TI = ndof*GlobalBasisFuncs[ii];

      b1 = SolnData->var1Cur(TI);
      b2 = SolnData->var1Cur(TI+1);
      b4 = SolnData->var1(TI+2);

      vel(0) += ( b1 * N(ii) );
      vel(1) += ( b2 * N(ii) );

      stress(0,0) += ( b1 * dN_dx(ii) );
      stress(0,1) += ( b1 * dN_dy(ii) );

      stress(1,0) += ( b2 * dN_dx(ii) );
      stress(1,1) += ( b2 * dN_dy(ii) );

      pres   += ( b4 * N(ii) );
    }

    // this is pseudo-stress
    stress *= mu;
    //stress = mu*(stress+stress.transpose());

    stress(0,0) -= pres;
    stress(1,1) -= pres;

  return;
}



template<>
void TreeNode<3>::computeVelocityAndStress(myPoint& param, myPoint& vel, MatrixXd& stress, myPoint& acce)
{
    int ii, TI, TIp1, TIp2, TIp3;

    double  b1, b2, b3, b4, pres;

    VectorXd  NN(totnlbf), dNN_dx(totnlbf), dNN_dy(totnlbf), dNN_dz(totnlbf);
    VectorXd  N, dN_dx, dN_dy, dN_dz;

    double  mu  = elmDat[4];

    GeomData->computeBasisFunctions3D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy, dNN_dz);

    if(parent == NULL)
    {
      N = NN;

      dN_dx = dNN_dx;
      dN_dy = dNN_dy;
      dN_dz = dNN_dz;
    }
    else
    {
      N = SubDivMat*NN;

      dN_dx = SubDivMat*dNN_dx;
      dN_dy = SubDivMat*dNN_dy;
      dN_dz = SubDivMat*dNN_dz;
    }

    vel.setZero();
    stress.setZero();
    pres = 0.0;

    for(ii=0;ii<totnlbf2;ii++)
    {
      TI = ndof*GlobalBasisFuncs[ii];

      b1 = SolnData->var1(TI);
      b2 = SolnData->var1(TI+1);
      b3 = SolnData->var1(TI+2);
      b4 = SolnData->var1(TI+3);

      vel(0) += ( b1 * N(ii) );
      vel(1) += ( b2 * N(ii) );
      vel(2) += ( b3 * N(ii) );

      stress(0,0) += ( b1 * dN_dx(ii) );
      stress(0,1) += ( b1 * dN_dy(ii) );
      stress(0,2) += ( b1 * dN_dz(ii) );

      stress(1,0) += ( b2 * dN_dx(ii) );
      stress(1,1) += ( b2 * dN_dy(ii) );
      stress(1,2) += ( b2 * dN_dz(ii) );

      stress(2,0) += ( b3 * dN_dx(ii) );
      stress(2,1) += ( b3 * dN_dy(ii) );
      stress(2,2) += ( b3 * dN_dz(ii) );

      pres   += ( b4 * N(ii) );
    }

    // this is pseudo-stress
    stress *= mu;
    //stress = mu*(stress+stress.transpose());

    stress(0,0) -= pres;
    stress(1,1) -= pres;
    stress(2,2) -= pres;

  return;
}


template<>
void TreeNode<3>::computeVelocityAndStressCur(myPoint& param, myPoint& vel, MatrixXd& stress, myPoint& acce)
{
    int ii, TI;

    double  b1, b2, b3, b4, pres;

    VectorXd  NN(totnlbf), dNN_dx(totnlbf), dNN_dy(totnlbf), dNN_dz(totnlbf);
    VectorXd  N, dN_dx, dN_dy, dN_dz;

    double  mu  = elmDat[4];

    GeomData->computeBasisFunctions3D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy, dNN_dz);

    if(parent == NULL)
    {
      N = NN;

      dN_dx = dNN_dx;
      dN_dy = dNN_dy;
      dN_dz = dNN_dz;
    }
    else
    {
      N = SubDivMat*NN;

      dN_dx = SubDivMat*dNN_dx;
      dN_dy = SubDivMat*dNN_dy;
      dN_dz = SubDivMat*dNN_dz;
    }

    vel.setZero();
    stress.setZero();
    pres = 0.0;

    for(ii=0;ii<totnlbf2;ii++)
    {
      TI = ndof*GlobalBasisFuncs[ii];

      b1 = SolnData->var1Cur(TI);
      b2 = SolnData->var1Cur(TI+1);
      b3 = SolnData->var1Cur(TI+2);
      b4 = SolnData->var1(TI+3);

      vel(0) += ( b1 * N(ii) );
      vel(1) += ( b2 * N(ii) );
      vel(2) += ( b3 * N(ii) );

      stress(0,0) += ( b1 * dN_dx(ii) );
      stress(0,1) += ( b1 * dN_dy(ii) );
      stress(0,2) += ( b1 * dN_dz(ii) );

      stress(1,0) += ( b2 * dN_dx(ii) );
      stress(1,1) += ( b2 * dN_dy(ii) );
      stress(1,2) += ( b2 * dN_dz(ii) );

      stress(2,0) += ( b3 * dN_dx(ii) );
      stress(2,1) += ( b3 * dN_dy(ii) );
      stress(2,2) += ( b3 * dN_dz(ii) );

      pres   += ( b4 * N(ii) );
    }

    // this is pseudo-stress
    stress *= mu;
    //stress = mu*(stress+stress.transpose());

    stress(0,0) -= pres;
    stress(1,1) -= pres;
    stress(2,2) -= pres;

  return;
}



