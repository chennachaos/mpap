
#include "TreeNode.h"
#include "MpapTime.h"
#include "Functions.h"
#include "TimeFunction.h"
#include "SolutionData.h"

#include "BasisFunctionsBSpline.h"

extern MpapTime mpapTime;
extern List<TimeFunction> timeFunction;


template<>
int TreeNode<2>::calcLoadVector(int ind1, int ind2, double inp1, double inp2)
{
   cout << " TreeNode<2>::calcLoadVector() ...... yet to be implemented " << endl;

  return 0;
}




/*
template<>
void TreeNode<2>::calcStiffnessAndResidualLSFEM(MatrixXd& Klocal, VectorXd& Flocal)
{
    // LSFEM for Stokes (transient) ---- new
    //
    //////////////////////////////////////////////////

    int ii, jj, gp1, gp2, TI, TIp1, TIp2, count, TJ, TJp1, TJp2;
   
    double  Da, Db, af, am, d1, c1, muTaf, rad, urdr, urdr2;
    double  uu, vv, Jac, dvol, fact, b1, b2, b3, fact1, fact2, fact3, acceFact;
    double  xx, yy, pres, HH, OMamdGamma, c, dist;

    af = SolnData->td(2);
    am = SolnData->td(3);
    acceFact = rho*am*SolnData->td(9);
    muTaf = mu*af;
    
    VectorXd  NN(totnlbf), dNN_dx(totnlbf), d2NN_dx2(totnlbf), dNN_dy(totnlbf), d2NN_dy2(totnlbf), vectmp(totnlbf);
    VectorXd  N, dN_dx, d2N_dx2, dN_dy, d2N_dy2;
    VectorXd  res(3), dp(2), Du(2), vel(2), velDot(2), force(2), R(2);
    MatrixXd  D(nsize2, 3), F(2,2), FN(2,2), stress(2,2);
    D.setZero();

    count = 0;
    for(gp2=0;gp2<GeomData->GetNGP(1);gp2++)
    {
       vv  = 0.5*(knots[1][2] * GeomData->gausspoints2[gp2] + knots[1][3]);
       Jac = GeomData->gaussweights2[gp2] * JacMult;
       
       for(gp1=0;gp1<GeomData->GetNGP(0);gp1++)
       {
          uu   = 0.5*(knots[0][2] * GeomData->gausspoints1[gp1] + knots[0][3]);
          dvol = GeomData->gaussweights1[gp1] * Jac;

          GeomData->computeBasisFunctions2D(knots[0][0], knots[1][0], knots[0][2], knots[1][2], uu, vv, &NN(0), &dNN_dx(0), &dNN_dy(0), &d2NN_dx2(0), &d2NN_dy2(0));

          //N = GeomData->shpfns[level][count].N;
          //dN_dx = GeomData->shpfns[level][count].dN_dx;
          //dN_dy = GeomData->shpfns[level][count].dN_dy;
          //d2N_dx2 = GeomData->shpfns[level][count].d2N_dx2;
          //d2N_dy2 = GeomData->shpfns[level][count].d2N_dy2;
          //count++;

          if(parent == NULL)
          {
            N = NN;
            dN_dx = dNN_dx;
            dN_dy = dNN_dy;
            d2N_dx2 = d2NN_dx2;
            d2N_dy2 = d2NN_dy2;
          }
          else
          {
            N = SubDivMat*NN;
            dN_dx = SubDivMat*dNN_dx;
            dN_dy = SubDivMat*dNN_dy;
            d2N_dx2 = SubDivMat*d2NN_dx2;
            d2N_dy2 = SubDivMat*d2NN_dy2;
          }

          //for(ii=0;ii<totnlbf;ii++)
            //printf(" \t %14.8f\t %14.8f\t%14.8f\t %14.8f\t%14.8f\n", N[ii], dN_dx[ii], dN_dy[ii], d2N_dx2[ii], d2N_dy2[ii]);

          //xx = GeomData->ComputeCoord(0, uu);
          //yy = GeomData->ComputeCoord(1, vv);

          //if(dist > 0.0)
            //HH = 1.0;
          //else
            //HH = 0.0;

          vel(0) = computeValueCur(0, N);
          vel(1) = computeValueCur(1, N);

          F(0,0) = computeValueCur(0, dN_dx);
          F(0,1) = computeValueCur(0, dN_dy);
          F(1,0) = computeValueCur(1, dN_dx);
          F(1,1) = computeValueCur(1, dN_dy);

          vectmp = d2N_dx2 + d2N_dy2;
          Du(0) = computeValueCur(0, vectmp);
          Du(1) = computeValueCur(1, vectmp);

          pres  = computeValue(2, N);
          dp(0) = computeValue(2, dN_dx);
          dp(1) = computeValue(2, dN_dy);

          velDot(0) = computeValueDotCur(0, N);
          velDot(1) = computeValueDotCur(1, N);

          for(ii=0;ii<totnlbf2;ii++)
          {
             TI   =  3*ii;
             TIp1 =  TI+1;
             TIp2 =  TI+2;

             b1 = dN_dx[ii];
             b2 = dN_dy[ii];
             b3 = rho*af*N[ii];

             fact = acceFact*N[ii] - mu*vectmp[ii];
             fact *= af;
             
             D(TI,0)   = fact ;
             D(TIp1,0) = 0.0;
             D(TIp2,0) = b1;

             D(TI,1)   = 0.0;
             D(TIp1,1) = fact;
             D(TIp2,1) = b2;
           
             D(TI,2)   = af*b1;
             D(TIp1,2) = af*b2;
             D(TIp2,2) = 0.0;
          }

          res.setZero();

          res(0) = computeForceCur(0, N) ;
          res(1) = computeForceCur(1, N) ;

          //res(0) = analy.computeForce(0, xx, yy);
          //res(1) = analy.computeForce(1, xx, yy);

          R = rho*velDot - mu*Du + dp;

          res(0) -= R(0);
          res(1) -= R(1);
          res(2) -= F.trace();

          //if(ind1)
          Klocal += ((dvol*D)*D.transpose());
          Flocal += (D*(dvol*res));
       }
    }
    //printf("\n\n");
    //printVector(Flocal);
    //printf("\n\n");
    //printMatrix(Klocal);
    //printf("\n\n");

   return;
}
*/


/*
template<>
void TreeNode<2>::computeAndReturnJacobian(int index, double* position, double* normal, double* specVal, double arclen, double* data, MatrixXd& D, VectorXd& vec1, VectorXd& vec2)
{
    // LSFEM for Stokes (transient)
    //
    //////////////////////////////////////////////////

    int ii, jj, gp1, gp2, TI, TIp1, TIp2, count, TJ;
    double  uu, vv, Jac, dvol, fact, fact1, fact2, fact3;
    double  pres, b1, b2, b3, xx, yy, afdt, af, HH, dt, am;

    dt = SolnData->td(0);
    af = SolnData->td(2);
    am = SolnData->td(3);
    //gamma = SolnData->td(4);
    afdt = SolnData->td(5);
    //c = SolnData->td(6);
    //OMamdGamma = SolnData->td(8);

    VectorXd  N(totnlbf), dN_dx(totnlbf), d2N_dx2(totnlbf), dN_dy(totnlbf), d2N_dy2(totnlbf), vectmp;
    VectorXd  res(2), dp(2), Du(2), vel(3), vel2(3);
    VectorXd  dpOld(2), DuOld(2), velOld(2), R(2), velCur(2), velDot(2);

    MatrixXd  D1(nsize, 2), D2(nsize, 2), F(2,2), NN(nsize,3), FOld(2,2), FCur(2,2);

    fact1 = 0.0;
    fact2 = 0.0;
    D2.setZero();
    NN.setZero();
    vel.setZero();

    count = 0;
    for(gp2=0;gp2<GeomData->GetNGP(1);gp2++)
    {
       vv  = 0.5*(knots[1][2] * GeomData->gausspoints2[gp2] + knots[1][3]);
       Jac = GeomData->gaussweights2[gp2] * JacMult;
       
       for(gp1=0;gp1<GeomData->GetNGP(0);gp1++)
       {
          uu   = 0.5*(knots[0][2] * GeomData->gausspoints1[gp1] + knots[0][3]);
          dvol = GeomData->gaussweights1[gp1] * Jac;

          //computeBasisFunctions2D(knots[0][0], knots[1][0], knots[0][2], knots[1][2], uu, vv, &N(0), &dN_dx(0), &dN_dy(0), &d2N_dx2(0), &d2N_dy2(0));

          N = GeomData->shpfns[level][count].N;
          dN_dx = GeomData->shpfns[level][count].dN_dx;
          dN_dy = GeomData->shpfns[level][count].dN_dy;
          d2N_dx2 = GeomData->shpfns[level][count].d2N_dx2;
          d2N_dy2 = GeomData->shpfns[level][count].d2N_dy2;
          count++;

          //for(ii=0;ii<totnlbf;ii++)
            //printf(" \t %14.8f\t %14.8f\t%14.8f\t %14.8f\t%14.8f\n", N[ii], dN_dx[ii], dN_dy[ii], d2N_dx2[ii], d2N_dy2[ii]);

          xx = GeomData->ComputeCoord(0, uu);
          yy = GeomData->ComputeCoord(1, vv);

          velCur(0) = computeValue(0, N);
          velCur(1) = computeValue(1, N);

          vectmp = d2N_dx2 + d2N_dy2;

          Du(0) = computeValue(0, vectmp);
          Du(1) = computeValue(1, vectmp);

          dp(0) = computeValue(2, dN_dx);
          dp(1) = computeValue(2, dN_dy);

          FCur(0,0) = computeValue(0, dN_dx);
          FCur(0,1) = computeValue(0, dN_dy);
          FCur(1,0) = computeValue(1, dN_dx);
          FCur(1,1) = computeValue(1, dN_dy);

          velOld(0) = computeValueOld(0, N);
          velOld(1) = computeValueOld(1, N);

          DuOld(0) = computeValueOld(0, vectmp);
          DuOld(1) = computeValueOld(1, vectmp);

          dpOld(0) = computeValueOld(2, dN_dx);
          dpOld(1) = computeValueOld(2, dN_dy);

          FOld(0,0) = computeValueOld(0, dN_dx);
          FOld(0,1) = computeValueOld(0, dN_dy);
          FOld(1,0) = computeValueOld(1, dN_dx);
          FOld(1,1) = computeValueOld(1, dN_dy);

          vel = af*velCur + (1.0-af)*velOld;
          F   = af*FCur   + (1.0-af)*FOld;

          for(ii=0;ii<totnlbf;ii++)
          {
             TI   =  3*ii;
             TIp1 =  TI+1;
             TIp2 =  TI+2;

             b1 = dN_dx[ii];
             b2 = dN_dy[ii];
             b3 = rho*afdt*N[ii];

             fact = rho*N[ii] - afdt*mu*vectmp[ii];
             
             D1(TI,0)   = fact; 
             D1(TIp1,0) = 0.0;
             D1(TIp2,0) = dt*b1;

             D1(TI,1)   = 0.0;
             D1(TIp1,1) = fact; 
             D1(TIp2,1) = dt*b2;
           
             NN(TI,0)   = N[ii];
             NN(TIp1,1) = N[ii];
          }

          //printf(" \t %14.8f\t %14.8f\t%14.8f\t %14.8f\t%14.8f\n", res(0), res(1), delta1, delta2, delta);

          delta *= arclen;

          fact = delta*dvol;

          D2 += ((fact*dt)*D1);

          //Ddelta += ((delta*fact)*(normal[0]*normal[0]+normal[1]*normal[1]));
          Ddelta += (delta*fact*dt*dt);

          velDot.setZero();
          res.setZero();

          //res(0) = computeForce(0, N) ;
          //res(1) = computeForce(1, N) ;
          //res *= dt;

          R = - mu*(af*Du+(1.0-af)*DuOld) + dp;

          R(0) += delta*data[0];
          R(1) += delta*data[1];

          res(0) = (velDot(0) + rho*(velCur(0)-velOld(0)) + dt*R(0) );
          res(1) = (velDot(1) + rho*(velCur(1)-velOld(1)) + dt*R(1) );

          fact1  += (0.0-res(0))*fact*dt; // 0.0 is the body force. need to generalise this
          fact2  += (0.0-res(1))*fact*dt;

          vel(0) = computeValue(0, N);
          vel(1) = computeValue(1, N);

          vel2.setZero();

          vel2(0) = 0.0 - vel(0);
          vel2(1) = 0.0 - vel(1);

          delta = delta1 * delta2;

          fact = delta*delta*dvol*1.0;
          Klocal += ( (fact*NN)*NN.transpose());
          Flocal += ( NN*(vel2*fact));
       }
    }

    //printf("\n\n");
    //printVector(Flocal);
    //printMatrix(D2);
    //printf("\n\n");

    if(parent == NULL)
    {
      D = D2;
    }
    else
    {
      D.resize(nsize2,2);
      D = SubDivMat2*D2;
    }

    vec1(0) = Ddelta;
    vec1(1) = Ddelta;
    vec2(0) = fact1;
    vec2(1) = fact2;

   return;
}
*/



template<>
void TreeNode<2>::calcResidualLSFEM(VectorXd& Flocal)
{
    // LSFEM for Navier-Stokes (transient) 3 DOF
    //
    //////////////////////////////////////////////////

    //Circle  circle2(10.0,15.0,0.5);

    int ii, jj, gp1, gp2, TI, TIp1, TIp2, count;
    double  Jac, dvol, fact, b1, b2, b3, fact1, fact2, fact3;
    double  xx, yy, pres, af, am, dist, acceFact, muTaf;

    count  = forAssyVec.size();

    bool   axsy = ((int)elmDat[2] == 1);
    double  rho = elmDat[3];
    double  mu  = elmDat[4];

    af = SolnData->td(2);
    am = SolnData->td(3);
    acceFact = am*SolnData->td(9);
    muTaf = mu*af;

    VectorXd  NN(totnlbf), dNN_dx(totnlbf), d2NN_dx2(totnlbf), dNN_dy(totnlbf), d2NN_dy2(totnlbf);
    VectorXd  N, dN_dx, d2N_dx2, dN_dy, d2N_dy2, vectmp;
    VectorXd  res(3), dp(2), Du(2), vel(2), velOld(2), R(2), velDot(2);
    MatrixXd  D(count, 3), F(2,2);
    myPoint  param;

    count = 0;
    for(gp2=0;gp2<GeomData->GetNGP(1);gp2++)
    {
       param[1]  = 0.5*(knots[1][2] * GeomData->gausspoints2[gp2] + knots[1][3]);
       Jac = GeomData->gaussweights2[gp2] * JacMultElem;
       
       for(gp1=0;gp1<GeomData->GetNGP(0);gp1++)
       {
          param[0]   = 0.5*(knots[0][2] * GeomData->gausspoints1[gp1] + knots[0][3]);
          dvol = GeomData->gaussweights1[gp1] * Jac;

          GeomData->computeBasisFunctions2D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy, d2NN_dx2, d2NN_dy2);

          if(parent == NULL)
          {
            N       = NN;
            dN_dx   = dNN_dx;
            dN_dy   = dNN_dy;
            d2N_dx2 = d2NN_dx2;
            d2N_dy2 = d2NN_dy2;
          }
          else
          {
            N       = SubDivMat*NN;
            dN_dx   = SubDivMat*dNN_dx;
            dN_dy   = SubDivMat*dNN_dy;
            d2N_dx2 = SubDivMat*d2NN_dx2;
            d2N_dy2 = SubDivMat*d2NN_dy2;
          }

          //for(ii=0;ii<totnlbf;ii++)
            //printf(" \t %14.8f\t %14.8f\t%14.8f\t %14.8f\t%14.8f\n", N[ii], dN_dx[ii], dN_dy[ii], d2N_dx2[ii], d2N_dy2[ii]);

          velDot(0) = computeValueDotCur(0, N);
          velDot(1) = computeValueDotCur(1, N);

          vel(0) = computeValueCur(0, N);
          vel(1) = computeValueCur(1, N);

          velOld(0) = computeValuePrev(0, N);
          velOld(1) = computeValuePrev(1, N);

          vectmp = d2N_dx2 + d2N_dy2;
          Du(0) = computeValueCur(0, vectmp);
          Du(1) = computeValueCur(1, vectmp);

          dp(0) = computeValue(2, dN_dx);
          dp(1) = computeValue(2, dN_dy);

          F(0,0) = computeValueCur(0, dN_dx);
          F(0,1) = computeValueCur(0, dN_dy);
          F(1,0) = computeValueCur(1, dN_dx);
          F(1,1) = computeValueCur(1, dN_dy);

          //printf(" \t %14.8f\t %14.8f\t%14.8f\t %14.8f\n", dist, HH, res(0), res(1));

          //xx = GeomData->ComputeCoord(0, uu);
          //yy = GeomData->ComputeCoord(1, vv);

          for(ii=0;ii<totnlbf2;ii++)
          {
             TI   =  3*ii;
             TIp1 =  TI+1;
             TIp2 =  TI+2;

             b1 = dN_dx[ii];
             b2 = dN_dy[ii];
             b3 = rho*af*N[ii];

             //fact = rho*acceFact*N[ii] + af*(rho*(velOld(0)*b1 + velOld(1)*b2) - mu*vectmp[ii]);
             fact = rho*acceFact*N[ii] + af*(rho*(vel(0)*b1 + vel(1)*b2) - mu*vectmp[ii]);

             D(TI,0)   = F(0,0)*b3 + fact ; 
             D(TIp1,0) = F(0,1)*b3;
             D(TIp2,0) = b1;

             D(TI,1)   = F(1,0)*b3;
             D(TIp1,1) = F(1,1)*b3 + fact;
             D(TIp2,1) = b2;

             D(TI,2)   = af*b1;
             D(TIp1,2) = af*b2;
             D(TIp2,2) = 0.0;
          }


          //R = rho*(velDot+F*velOld) - mu*Du + dp;
          R = rho*(velDot+F*vel) - mu*Du + dp;

          res.setZero();
          res(0) = computeForceCur(0, N) ;
          res(1) = computeForceCur(1, N) ;

          res(0) -= R(0) ;
          res(1) -= R(1) ;
          res(2) -= F.trace() ;

          Flocal += (D*(dvol*res));
       }
    }
    
   return;
}



//
template<>
void TreeNode<2>::calcStiffnessAndResidualLSFEM(bool flag, MatrixXd& Klocal, VectorXd& Flocal)
{
    // LSFEM for Navier-Stokes (transient) 4 DOF
    //
    //////////////////////////////////////////////////

    int ii, jj, gp1, gp2, TI, TIp1, TIp2, TIp3, count;
    double  uu, vv, Jac, dvol, fact, b1, b2, b3, fact1, fact2, fact3;
    double  xx, yy, pres, af, am, HH, dist, w, acceFact, muTaf, dt;

    bool   axsy = ((int)elmDat[2] == 1);
    double  rho = elmDat[3];
    double  mu  = elmDat[4];


    dt = SolnData->td(0);
    af = SolnData->td(2);
    am = SolnData->td(3);
    acceFact = am*SolnData->td(9);

    count  = forAssyVec.size();

    //cout << totnlbf << '\t' << totnlbf2 << '\t' << nU << '\t' << nP << '\t' << count << endl;

    VectorXd  NN(totnlbf), dNN_dx(totnlbf), dNN_dy(totnlbf);
    VectorXd  N, dN_dx, dN_dy;
    VectorXd  res(ndof), dp(2), vel(2), dw(2), R(2), velDot(2), velOld(2);
    MatrixXd  D(count, ndof), F(2,2);
    myPoint  param;

    count = 0;
    for(gp2=0;gp2<GeomData->GetNGP(1);gp2++)
    {
       param[1]  = 0.5*(knots[1][2] * GeomData->gausspoints2[gp2] + knots[1][3]);
       Jac = GeomData->gaussweights2[gp2] * JacMultElem;
       
       for(gp1=0;gp1<GeomData->GetNGP(0);gp1++)
       {
          param[0]  = 0.5*(knots[0][2] * GeomData->gausspoints1[gp1] + knots[0][3]);
          dvol = GeomData->gaussweights1[gp1] * Jac;

          GeomData->computeBasisFunctions2D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy);

          if(parent == NULL)
          {
            N       = NN;
            dN_dx   = dNN_dx;
            dN_dy   = dNN_dy;
          }
          else
          {
            N       = SubDivMat*NN;
            dN_dx   = SubDivMat*dNN_dx;
            dN_dy   = SubDivMat*dNN_dy;
          }

          //for(ii=0;ii<totnlbf;ii++)
            //printf(" \t %14.8f\t %14.8f\t%14.8f\t %14.8f\t%14.8f\n", N[ii], dN_dx[ii], dN_dy[ii], d2N_dx2[ii], d2N_dy2[ii]);

          //xx = GeomData->ComputeCoord(0, uu);
          //yy = GeomData->ComputeCoord(1, vv);

          //pt[0] = xx;
          //pt[1] = yy;

          velDot(0) = computeValueDotCur(0, N);
          velDot(1) = computeValueDotCur(1, N);

          vel(0) = computeValueCur(0, N);
          vel(1) = computeValueCur(1, N);

          velOld(0) = computeValuePrev(0, N);
          velOld(1) = computeValuePrev(1, N);

          F(0,0) = computeValueCur(0, dN_dx);
          F(0,1) = computeValueCur(0, dN_dy);
          F(1,0) = computeValueCur(1, dN_dx);
          F(1,1) = computeValueCur(1, dN_dy);

          dp(0) = computeValue(2, dN_dx);
          dp(1) = computeValue(2, dN_dy);

          w      = computeValueCur(3, N);
          dw(0)  = computeValueCur(3, dN_dx);
          dw(1)  = computeValueCur(3, dN_dy);

          //printf(" \t %14.8f\t %14.8f\t%14.8f\t %14.8f\n", dist, HH, res(0), res(1));

          for(ii=0;ii<totnlbf2;ii++)
          {
             TI   =  ndof*ii;
             TIp1 =  TI+1;
             TIp2 =  TI+2;
             TIp3 =  TI+3;

             b1 = dN_dx[ii];
             b2 = dN_dy[ii];
             b3 = rho*af*dt*N[ii];

             fact = rho*dt*acceFact*N[ii] + rho*af*dt*(vel(0)*b1 + vel(1)*b2);

             D(TI,0)   = F(0,0)*b3 + fact;
             D(TIp1,0) = F(0,1)*b3 ;
             D(TIp2,0) = dt*b1;
             D(TIp3,0) = mu*af*dt*b2;

             D(TI,1)   = F(1,0)*b3 ;
             D(TIp1,1) = F(1,1)*b3 + fact;
             D(TIp2,1) = dt*b2;
             D(TIp3,1) = -mu*af*dt*b1;

             //fact = rho*dt*acceFact*N[ii] + rho*af*dt*(velOld(0)*b1 + velOld(1)*b2);

             //D(TI,0)   = fact;
             //D(TIp1,0) = 0.0 ;
             //D(TIp2,0) = dt*b1;
             //D(TIp3,0) = mu*af*dt*b2;

             //D(TI,1)   = 0.0 ;
             //D(TIp1,1) = fact;
             //D(TIp2,1) = dt*b2;
             //D(TIp3,1) = -mu*af*dt*b1;

             D(TI,2)   =  af*b1;
             D(TIp1,2) =  af*b2;
             D(TIp2,2) =  0.0;
             D(TIp3,2) =  0.0;

             D(TI,3)   =  af*b2;
             D(TIp1,3) = -af*b1;
             D(TIp2,3) = 0.0;
             D(TIp3,3) = af*N(ii);
          }

          res.setZero();

          //res(0) = dt*computeForceCur(0, N) ;
          //res(1) = dt*computeForceCur(1, N) ;

          R = rho*(velDot+F*vel) + dp;

          res(0) -= dt*( R(0) + mu*dw(1) );
          res(1) -= dt*( R(1) - mu*dw(0) );
          res(2) -= (F.trace());
          res(3) -= (w+F(0,1)-F(1,0));

          //printf("\n\n \t %14.8f\t %14.8f\t%14.8f\t %14.8f\t%14.8f\t%14.8f\t%14.8f\t%14.8f\t%14.8f\n\n\n", res(0), res(1), F(0,0), F(0,1), F(1,0), F(1,1), JacMult, Jac, dvol);
          //printf(" \t %14.8f\t %14.8f\t%14.8f\t %14.8f\t%14.8f\t%14.8f\t%14.8f\n", res(0), res(1), res(2), F(0,0), F(0,1), F(1,0), F(1,1));

          //if(ind1)
          Klocal += ((dvol*D)*D.transpose());
          Flocal += (D*(dvol*res));
       }
    }

   return;
}
//

          /*
          if(parent == NULL)
          {
            N       = GeomData->shpfns[level][count].N;
            dN_dx   = GeomData->shpfns[level][count].dN_dx;
            dN_dy   = GeomData->shpfns[level][count].dN_dy;
            d2N_dx2 = GeomData->shpfns[level][count].d2N_dx2;
            d2N_dy2 = GeomData->shpfns[level][count].d2N_dy2;
          }
          else
          {
            N       = SubDivMat*GeomData->shpfns[level][count].N;
            dN_dx   = SubDivMat*GeomData->shpfns[level][count].dN_dx;
            dN_dy   = SubDivMat*GeomData->shpfns[level][count].dN_dy;
            d2N_dx2 = SubDivMat*GeomData->shpfns[level][count].d2N_dx2;
            d2N_dy2 = SubDivMat*GeomData->shpfns[level][count].d2N_dy2;
          }
          count++;
          */


/*
template<>
void TreeNode<2>::calcStiffnessAndResidualLSFEM(bool flag, MatrixXd& Klocal, VectorXd& Flocal)
{
    // LSFEM for Navier-Stokes (transient) 3 DOF
    //
    //////////////////////////////////////////////////

    //Circle  circle2(10.0,15.0,0.5);

    int ii, jj, gp1, gp2, TI, TIp1, TIp2, count;
    double  uu, vv, Jac, dvol, fact, b1, b2, b3, fact1, fact2, fact3;
    double  xx, yy, pres, af, am, dist, acceFact, muTaf, dt;

    count  = forAssyVec.size();

    dt = SolnData->td(0);
    af = SolnData->td(2);
    am = SolnData->td(3);
    acceFact = am*SolnData->td(9);

    VectorXd  NN(totnlbf), dNN_dx(totnlbf), d2NN_dx2(totnlbf), dNN_dy(totnlbf), d2NN_dy2(totnlbf);
    VectorXd  N, dN_dx, d2N_dx2, dN_dy, d2N_dy2, vectmp;
    VectorXd  res(3), dp(2), Du(2), vel(2), velOld(2), R(2), velDot(2);
    MatrixXd  D(count, 3), F(2,2);

    count = 0;
    for(gp2=0;gp2<GeomData->GetNGP(1);gp2++)
    {
       vv  = 0.5*(knots[1][2] * GeomData->gausspoints2[gp2] + knots[1][3]);
       Jac = GeomData->gaussweights2[gp2] * JacMult;
       
       for(gp1=0;gp1<GeomData->GetNGP(0);gp1++)
       {
          uu   = 0.5*(knots[0][2] * GeomData->gausspoints1[gp1] + knots[0][3]);
          dvol = GeomData->gaussweights1[gp1] * Jac;

          GeomData->computeBasisFunctions2D(knots[0][0], knots[1][0], knots[0][2], knots[1][2], uu, vv, &NN(0), &dNN_dx(0), &dNN_dy(0), &d2NN_dx2(0), &d2NN_dy2(0));

          if(parent == NULL)
          {
            N       = NN;
            dN_dx   = dNN_dx;
            dN_dy   = dNN_dy;
            d2N_dx2 = d2NN_dx2;
            d2N_dy2 = d2NN_dy2;
          }
          else
          {
            N       = SubDivMat*NN;
            dN_dx   = SubDivMat*dNN_dx;
            dN_dy   = SubDivMat*dNN_dy;
            d2N_dx2 = SubDivMat*d2NN_dx2;
            d2N_dy2 = SubDivMat*d2NN_dy2;
          }

          //for(ii=0;ii<totnlbf;ii++)
            //printf(" \t %14.8f\t %14.8f\t%14.8f\t %14.8f\t%14.8f\n", N[ii], dN_dx[ii], dN_dy[ii], d2N_dx2[ii], d2N_dy2[ii]);

          velDot(0) = computeValueDotCur(0, N);
          velDot(1) = computeValueDotCur(1, N);

          vel(0) = computeValueCur(0, N);
          vel(1) = computeValueCur(1, N);

          velOld(0) = computeValueOld(0, N);
          velOld(1) = computeValueOld(1, N);

          vectmp = d2N_dx2 + d2N_dy2;
          Du(0) = computeValueCur(0, vectmp);
          Du(1) = computeValueCur(1, vectmp);

          dp(0) = computeValue(2, dN_dx);
          dp(1) = computeValue(2, dN_dy);

          F(0,0) = computeValueCur(0, dN_dx);
          F(0,1) = computeValueCur(0, dN_dy);
          F(1,0) = computeValueCur(1, dN_dx);
          F(1,1) = computeValueCur(1, dN_dy);

          //printf(" \t %14.8f\t %14.8f\t%14.8f\t %14.8f\n", dist, HH, res(0), res(1));

          for(ii=0;ii<totnlbf2;ii++)
          {
             TI   =  3*ii;
             TIp1 =  TI+1;
             TIp2 =  TI+2;

             b1 = dN_dx[ii];
             b2 = dN_dy[ii];
             b3 = rho*af*N[ii];

             //fact = rho*acceFact*N[ii] + af*(rho*(vel(0)*b1 + vel(1)*b2) - mu*vectmp[ii]);

             //D(TI,0)   = F(0,0)*b3 + fact ; 
             //D(TIp1,0) = F(0,1)*b3;
             //D(TIp2,0) = b1;

             //D(TI,1)   = F(1,0)*b3;
             //D(TIp1,1) = F(1,1)*b3 + fact;
             //D(TIp2,1) = b2;

             fact = rho*dt*acceFact*N[ii] + af*dt*(rho*(velOld(0)*b1 + velOld(1)*b2) - mu*vectmp[ii]);

             D(TI,0)   = fact ; 
             D(TIp1,0) = 0.0;
             D(TIp2,0) = dt*b1;

             D(TI,1)   = 0.0;
             D(TIp1,1) = fact;
             D(TIp2,1) = dt*b2;

             D(TI,2)   = rho*af*b1;
             D(TIp1,2) = rho*af*b2;
             D(TIp2,2) = 0.0;
          }


          R = rho*(velDot+F*velOld) - mu*Du + dp;
          //R = rho*(velDot+F*vel) - mu*Du + dp;

          res.setZero();
          res(0) = dt*computeForceCur(0, N) ;
          res(1) = dt*computeForceCur(1, N) ;

          res(0) -= dt*R(0) ;
          res(1) -= dt*R(1) ;
          res(2) -= rho*F.trace() ;

          //printf("\n\n \t %14.8f\t %14.8f\t%14.8f\t %14.8f\t%14.8f\t%14.8f\t%14.8f\t%14.8f\t%14.8f\n\n\n", res(0), res(1), F(0,0), F(0,1), F(1,0), F(1,1), JacMult, Jac, dvol);
          //printf(" \t %14.8f\t %14.8f\t%14.8f\t %14.8f\t%14.8f\t%14.8f\t%14.8f\n", res(0), res(1), res(2), F(0,0), F(0,1), F(1,0), F(1,1));

          //if(flag)
            Klocal += ((dvol*D)*D.transpose());

          Flocal += (D*(dvol*res));
       }
    }
    //printf("\n\n");
    //printVector(Flocal);
    //printf("\n\n");
    //printMatrix(Klocal);
    //printf("\n\n");
    
   return;
}
*/


template<>
void TreeNode<2>::applyDirichletBCsLSFEM(bool flag, MatrixXd& Klocal, VectorXd& Flocal)
{
  if( DirichletData.size() > 0 )
  {
    int ii, jj, aa, gp1, gp2, TI, TIp1, TIp2, index, TJ, TJp1, TJp2, dir, side, nU, nP;
    double theta, y0, y1, mu1, sigma, Ta, Tb, res, JacMultLoc, af, temp;
    double  dvol, specVal, PENALTY, xx, yy, r, Jac, fact, fact1, fact2;

    myPoint param, normal, trac;
    VectorXd  N(totnlbf), dN_dx(totnlbf), dN_dy(totnlbf);
    VectorXd  NN(totnlbf), dNN_dx(totnlbf), dNN_dy(totnlbf);
    vector<double>  boundaryGPs1, boundaryGWs1, boundaryGPs2, boundaryGWs2;

    bool   axsy = ((int)elmDat[2] == 1);
    double  rho = elmDat[3];
    double  mu  = elmDat[4];


    af = SolnData->td(2);
    
    y0 = 0.25;
    y1 = 0.75;
    theta = 0.0;
    //theta = 0.463647609;
    //y0 = 0.0;
    //y1 = 0.5;

    //y0 = 0.1;
    //y1 = 0.6;
    //theta = 0.291456794;

    //y0 = 0.0;
    //y1 = 1.0;

    //y0 = 0.0;
    //y1 = 1.61;

    //Kovasznay  analy;
    //analy.SetPressure(1.310741966654558776639305506250821053981781005859375);
    //analy.SetPressure(0.0);

    Stokes2DEx1  analy;
    //Stokes2DEx2  analy;

    //TwoDim_Ex1  analy;
    //PoissonEx2 analy;

    //PearsonVortex analy;

      for(aa=0;aa<DirichletData.size();aa++)
      {
        //printVector(DirichletData[aa]);

        side    = (int) (DirichletData[aa][0] - 1);
        dir     = (int) (DirichletData[aa][1] - 1);
        specVal = DirichletData[aa][2];
        PENALTY = DirichletData[aa][3];
        
        GeomData->getBoundaryNormal2D(side, normal);
        GeomData->setBoundaryGPs2D(side, boundaryGPs1, boundaryGWs1, boundaryGPs2, boundaryGWs2);

        JacMultLoc = TreeNode<2>::getJacBoundary(side);

        for(gp2=0;gp2<boundaryGPs2.size();gp2++)
        {
            param[1] = 0.5 * (knots[1][2] * boundaryGPs2[gp2] + knots[1][3]);
        for(gp1=0;gp1<boundaryGPs1.size();gp1++)
        {
            param[0] = 0.5 * (knots[0][2] * boundaryGPs1[gp1] + knots[0][3]);
 
            dvol = JacMultLoc * boundaryGWs2[gp2] * boundaryGWs1[gp1] ;
            
            //printf(" %4d \t %4d \t %12.6f \t %12.6f \n", side, dir, vv, uu);

            GeomData->computeBasisFunctions2D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy );

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

            //for(ii=0;ii<totnlbf;ii++)
            //printf(" \t %14.8f \n", NN[ii]);

            //printf("\n\n tracX and tracY ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n\n\n", xx, yy, val, val, Jac, JacMult, dvol);
            //if(pout) printf(" knotsAtGPs = %12.6f \t xx = %12.6f \t yy = %12.6f \n", knotsAtGPs, xx, yy);
            //printf(" uu = %12.6f \t vv = %12.6f \t dvol = %12.6f \t volume = %12.6f \n", uu, vv, dvol, volume);

            //xx = GeomData->ComputeCoord(0, uu);
            //yy = GeomData->ComputeCoord(1, vv);
            //r = sqrt(xx*xx+yy*yy);
            //val = 1.0+log(2.0*r);
            //cout << xx << '\t' << yy << endl;

            //cout << aa << '\t' << side << '\t' << dir << '\t' << val << endl;
            specVal = DirichletData[aa][2];

            /*
            if(side == 0)
            {
              if(dir == 0)
              {
                if(yy < y0 || yy > y1)
                  specVal = 0.0;
                else
                  specVal = DirichletData[aa][2]*(y1-yy)*(yy-y0);
              }
            }
            */

            //specVal = analy.computeSolution(dir, xx, yy);

            specVal *= timeFunction[0].prop;
            //specVal *= mpapTime.cur;
            //specVal *= tanh(2.0*mpapTime.cur);
            //specVal *= sin(2.0*PI*mpapTime.cur+PI/2.0);
            //specVal *= sin(2.0*PI*10.0*mpapTime.cur+PI/2.0);

            //cout << aa << '\t' << side << '\t' << dir << '\t' << res << endl;
            //cout << specVal << '\t' << computeValue(dir, N) << endl;

            res = specVal - computeValue(dir, N);

            for(ii=0;ii<totnlbf2;ii++)
            {
              fact = N[ii] * dvol * PENALTY;

              TI = ndof*ii+dir;

              Flocal(TI) += fact*res;
              //fact *= af;

              //if(flag)
              //{
                for(jj=0;jj<totnlbf2;jj++)
                  Klocal(TI, ndof*jj+dir) += fact * N[jj];
              //}
            }
        }// for(gp1=0...
        }// for(gp2=0...
      } // for(aa=0;aa<DirichletData.size();aa++)
  } // if(DirichletData.size() > 0)

  return;
}



/*
template<>
void TreeNode<2>::applyNeumannBCsLSFEM(bool flag, MatrixXd& Klocal, VectorXd& Flocal)
{
  if( NeumannData.size() > 0 )
  {
      int ii, aa, gp1, gp2, TI, TIp1, TIp2, dir, side, count;
      double  res, JacMult, af, uu, vv, dvol, specVal, PENALTY, xx, yy, pres;

      count  = forAssyVec.size();

      Point normal(2);
      VectorXd  N(totnlbf), dN_dx(totnlbf), dN_dy(totnlbf);
      VectorXd  NN(totnlbf), dNN_dx(totnlbf), dNN_dy(totnlbf);
      MatrixXd  D(count,1), stress(2,2);
      D.setZero();

      af = SolnData->td(2);

      PENALTY  = 1.0;

      for(aa=0;aa<NeumannData.size();aa++)
      {
        side    = (int) (NeumannData[aa][0] - 1);
        dir     = (int) (NeumannData[aa][1] - 1);
        specVal = NeumannData[aa][2];

        GeomData->setBoundaryGPs2D(side);
        GeomData->getBoundaryNormal2D(side, normal);

        JacMult = TreeNode<2>::getJacBoundary(side);

        for(gp2=0;gp2<GeomData->boundaryGPs2.size();gp2++)
        {
            vv = 0.5 * (knots[1][2] * GeomData->boundaryGPs2[gp2] + knots[1][3]);
        for(gp1=0;gp1<GeomData->boundaryGPs1.size();gp1++)
        {
            uu = 0.5 * (knots[0][2] * GeomData->boundaryGPs1[gp1] + knots[0][3]);
 
            dvol = JacMult * GeomData->boundaryGWs2[gp2] * GeomData->boundaryGWs1[gp1] ;

            GeomData->computeBasisFunctions2D(knots[0][0], knots[1][0], knots[0][2], knots[1][2], uu, vv, &NN(0), &dNN_dx(0), &dNN_dy(0) );

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

            //for(ii=0;ii<totnlbf;ii++)
            //printf(" \t %14.8f \n", NN[ii]);

            //printf("\n\n tracX and tracY ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n\n\n", xx, yy, val, val, Jac, JacMult, dvol);
            //if(pout) printf(" knotsAtGPs = %12.6f \t xx = %12.6f \t yy = %12.6f \n", knotsAtGPs, xx, yy);
            //printf(" uu = %12.6f \t vv = %12.6f \t dvol = %12.6f \t volume = %12.6f \n", uu, vv, dvol, volume);

            xx = GeomData->ComputeCoord(0, uu);
            yy = GeomData->ComputeCoord(1, vv);
            //r = sqrt(xx*xx+yy*yy);
            //val = 1.0+log(2.0*r);
            //cout << xx << '\t' << yy << endl;

            //specVal *= timeFunction[0].prop;
            //cout << val << '\t' << computeValue(dir, N) << endl;

            //res *= timeFunction[0].prop;
            //res *= (0.5*(1.0-cos(2.0*PI*SolnData->ElemProp.data[6]*mpapTime.cur)));
            //w1 = tanh(SolnData->ElemProp.data[5]*mpapTime.cur);
            //res *= w1;
            //w1 = 1.0;
            //w2 = 2.0*PI*SolnData->ElemProp.data[6];
            //res *= sin(w2*mpapTime.cur);
            //res *= tanh(50.0*timeFunction[0].prop);

            //stress(0,0) = 2.0*mu*computeValue(0, dN_dx);
            //stress(0,1) = mu*(computeValue(0, dN_dy)+computeValue(1, dN_dx));
            //stress(1,0) = stress(0,1);
            //stress(1,1) = 2.0*mu*computeValue(1, dN_dy);

            stress(0,0) = 2.0*mu*computeValueCur(0, dN_dx);
            stress(0,1) = mu*(computeValueCur(0, dN_dy)+computeValueCur(1, dN_dx));
            stress(1,0) = stress(0,1);
            stress(1,1) = 2.0*mu*computeValueCur(1, dN_dy);

            pres   = computeValue(2, N);

            stress(0,0) -= pres;
            stress(1,1) -= pres;

            specVal = NeumannData[aa][2];

            if(dir == 0)
            {
              for(ii=0;ii<totnlbf;ii++)
              {
                TI = ndof*ii;
                D(TI)   = 2.0*mu*af*dN_dx[ii]*normal[0] + mu*af*dN_dy[ii]*normal[1];
                D(TI+1) = mu*af*dN_dx[ii]*normal[1];
                D(TI+2) = -N[ii]*normal[0];
              }

              res = specVal - (stress(0,0)*normal[0]+stress(0,1)*normal[1]);
            }
            if(dir == 1)
            {
              for(ii=0;ii<totnlbf;ii++)
              {
                TI = ndof*ii;
                D(TI)   = mu*af*dN_dy[ii]*normal[0];
                D(TI+1) = mu*af*dN_dx[ii]*normal[0] + 2.0*mu*af*dN_dy[ii]*normal[1];
                D(TI+2) = -N[ii]*normal[1];
              }

              res = specVal - (stress(1,0)*normal[0]+stress(1,1)*normal[1]);
            }

            res *= dvol;

            Flocal += (res)*D;
            if(flag)
              Klocal += (dvol*D)*D.transpose();
        }// for(gp1=0...
        }// for(gp2=0...
      } // for(aa=0;aa<NeumannData.size();aa++)
  } // if(NeumannData.size() > 0)

  return;
}
*/



//
template<>
void TreeNode<2>::applyNeumannBCsLSFEM(bool flag, MatrixXd& Klocal, VectorXd& Flocal)
{
  if( NeumannData.size() > 0 )
  {
      int ii, aa, gp1, gp2, TI, TIp1, TIp2, dir, side, count;
      double  res, JacMultLoc, af, uu, vv, dvol, specVal, PENALTY, xx, yy, pres;
      myPoint  param;
      vector<double>   boundaryGPs1, boundaryGWs1, boundaryGPs2, boundaryGWs2;

      bool   axsy = ((int)elmDat[2] == 1);
      double  rho = elmDat[3];
      double  mu  = elmDat[4];

      count  = forAssyVec.size();
      cout << " count = " << count << endl;

      myPoint normal;
      VectorXd  N(totnlbf), dN_dx(totnlbf), dN_dy(totnlbf);
      VectorXd  NN(totnlbf), dNN_dx(totnlbf), dNN_dy(totnlbf);
      MatrixXd  D(count,1), stress(2,2);
      D.setZero();

      af = SolnData->td(2);

      PENALTY  = 1.0;

      for(aa=0;aa<NeumannData.size();aa++)
      {
        side    = (int) (NeumannData[aa][0] - 1);
        dir     = (int) (NeumannData[aa][1] - 1);
        specVal = NeumannData[aa][2];

        GeomData->setBoundaryGPs2D(side, boundaryGPs1, boundaryGWs1, boundaryGPs2, boundaryGWs2);
        GeomData->getBoundaryNormal2D(side, normal);

        JacMultLoc = TreeNode<2>::getJacBoundary(side);

        for(gp2=0;gp2<boundaryGPs2.size();gp2++)
        {
            param[1] = 0.5 * (knots[1][2] * boundaryGPs2[gp2] + knots[1][3]);
        for(gp1=0;gp1<boundaryGPs1.size();gp1++)
        {
            param[0] = 0.5 * (knots[0][2] * boundaryGPs1[gp1] + knots[0][3]);
 
            dvol = JacMultLoc * boundaryGWs2[gp2] * boundaryGWs1[gp1] ;

            GeomData->computeBasisFunctions2D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy );

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

            //for(ii=0;ii<totnlbf;ii++)
            //printf(" \t %14.8f \n", NN[ii]);

            //printf("\n\n tracX and tracY ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n\n\n", xx, yy, val, val, Jac, JacMult, dvol);
            //if(pout) printf(" knotsAtGPs = %12.6f \t xx = %12.6f \t yy = %12.6f \n", knotsAtGPs, xx, yy);
            //printf(" uu = %12.6f \t vv = %12.6f \t dvol = %12.6f \t volume = %12.6f \n", uu, vv, dvol, volume);

            xx = GeomData->ComputeCoord(0, uu);
            yy = GeomData->ComputeCoord(1, vv);
            //r = sqrt(xx*xx+yy*yy);
            //val = 1.0+log(2.0*r);
            cout << xx << '\t' << yy << endl;

            //specVal *= timeFunction[0].prop;
            //cout << val << '\t' << computeValue(dir, N) << endl;

            //res *= timeFunction[0].prop;
            //res *= (0.5*(1.0-cos(2.0*PI*SolnData->ElemProp.data[6]*mpapTime.cur)));
            //w1 = tanh(SolnData->ElemProp.data[5]*mpapTime.cur);
            //res *= w1;
            //w1 = 1.0;
            //w2 = 2.0*PI*SolnData->ElemProp.data[6];
            //res *= sin(w2*mpapTime.cur);
            //res *= tanh(50.0*timeFunction[0].prop);
            
            //stress(0,0) = computeValueCur(0, dN_dx);
            //stress(0,1) = computeValueCur(0, dN_dy);
            //stress(1,0) = computeValueCur(1, dN_dx);
            //stress(1,1) = computeValueCur(1, dN_dy);

            stress(0,0) = computeValue(0, dN_dx);
            stress(0,1) = computeValue(0, dN_dy);
            stress(1,0) = computeValue(1, dN_dx);
            stress(1,1) = computeValue(1, dN_dy);

            pres   = computeValue(2, N);
            stress = mu*stress;
            stress(0,0) -= pres;
            stress(1,1) -= pres;

            specVal = NeumannData[aa][2];
            //specVal = 0.0;

            if(dir == 0)
            {
              for(ii=0;ii<totnlbf2;ii++)
              {
                TI = ndof*ii;
                //D(TI)   = mu*af*(dN_dx[ii]*normal[0] + dN_dy[ii]*normal[1]);
                D(TI)   = mu*(dN_dx[ii]*normal[0] + dN_dy[ii]*normal[1]);
                D(TI+1) = 0.0;
                D(TI+2) = -N[ii]*normal[0];
              }

              res = specVal - (stress(0,0)*normal[0]+stress(0,1)*normal[1]);
            }
            if(dir == 1)
            {
              for(ii=0;ii<totnlbf2;ii++)
              {
                TI = ndof*ii;
                D(TI)   = 0.0;
                //D(TI+1) = mu*af*(dN_dx[ii]*normal[0] + dN_dy[ii]*normal[1]);
                D(TI+1) = mu*(dN_dx[ii]*normal[0] + dN_dy[ii]*normal[1]);
                D(TI+2) = -N[ii]*normal[1];
              }

              res = specVal - (stress(1,0)*normal[0]+stress(1,1)*normal[1]);
            }

            res *= dvol;

            Flocal += (res)*D;
            //if(flag)
              Klocal += (dvol*D)*D.transpose();
        }// for(gp1=0...
        }// for(gp2=0...
      } // for(aa=0;aa<NeumannData.size();aa++)
  } // if(NeumannData.size() > 0)
  cout << " AAAAAAAAAAA " << endl;

  return;
}
//




template<>
void TreeNode<2>::mapBoundaryPointDataToGlobalBodyForceVector(double* position, double* normal, double arclen, double* useThisData)
{
    int      ii, jj, gp1, gp2, TI;
    double   Jac, dvol, fact, xx, yy;
    double   beta1, beta2, delta1, delta2, delta, h1, h2;

    VectorXd  N(totnlbf), res(ndof);
    MatrixXd  D(nsize, ndof);
    myPoint  param;

    double  val1  = 0.5*knots[0][2];
    double  val2  = 0.5*knots[0][3];
    double  val3  = 0.5*knots[1][2];
    double  val4  = 0.5*knots[1][3];
    double  *gausspoints1  = &(GeomData->gausspoints1[0]);
    double  *gaussweights1 = &(GeomData->gaussweights1[0]);
    double  *gausspoints2  = &(GeomData->gausspoints2[0]);
    double  *gaussweights2 = &(GeomData->gaussweights2[0]);

    //JacMult = GeomData->GetJacobianFull() * val1 * val3;

    beta1 = elmDat[5] * GeomData->GetGridLength(0) * knots[0][2];
    beta2 = elmDat[5] * GeomData->GetGridLength(1) * knots[1][2];

    //h1 = GeomData->ElemProp.data[5] * GeomData->GetGridLength(0)*knots[0][2];
    //h2 = GeomData->ElemProp.data[5] * GeomData->GetGridLength(1)*knots[1][2];

    h1 = GeomData->GetGridLength(0) * knots[0][2];
    h2 = GeomData->GetGridLength(1) * knots[1][2];

    //printf("\t data \t %12.6f \t %12.6f \n", h1, h2);
    
    fact = arclen/h1/h2;
    fact = 1.0/h1/h2;
    //fact = 1.0;
    
    res.setZero();
    
    if(ndof == 1)
      res(0) = useThisData[0];
    else
    {
      res(0) = useThisData[0];
      res(1) = useThisData[1];
    }

    for(gp2=0;gp2<GeomData->GetNGP(1);gp2++)
    {
       param[1]  = val3 * gausspoints2[gp2] + val4;
       Jac = gaussweights2[gp2] * JacMultElem;

       for(gp1=0;gp1<GeomData->GetNGP(0);gp1++)
       {
          param[0]  = val1 * gausspoints1[gp1] + val2;
          dvol = gaussweights1[gp1] * Jac;

          GeomData->computeBasisFunctions2D(knotBegin, knotIncr, param, N);

          //printf(" \t %14.8f \t %14.8f \t %14.8f \t %14.8f\n", uu, vv, fact, dvol0);
          //printf("BasisFuns \n");
          //for(ii=0;ii<totnlbf;ii++)
            //printf(" \t %12.6f  \t %12.6f  \t %12.6f \n ", N(ii), dN_dx(ii), dN_dy(ii));

          xx = GeomData->ComputeCoord(0, param[0]) - position[0];
          yy = GeomData->ComputeCoord(1, param[1]) - position[1];

          delta1 = DiracDelta1(xx, beta1);
          delta2 = DiracDelta1(yy, beta2);

          delta = delta1 * delta2;

          dvol *= (delta*fact);

          for(ii=0;ii<totnlbf;ii++)
          {
            TI = ndof*ii;

            for(jj=0;jj<ndof;jj++)
              D(TI+jj,jj) = N[ii];
          }

          //if(delta > 0.0)
            //printf(" \t %14.8f\t %14.8f\t%14.8f\t %14.8f\t%14.8f\n", xx, yy, delta1, delta2, delta);

          //val1 = LevelSetFunc[0][7]*normal[0]*fact;
          //val2 = LevelSetFunc[0][7]*normal[1]*fact;

          Flocal += (dvol*(D*res));

    }//gp1
    }//gp2
    
    //printVector(Flocal);
    //printf("\n\n");

  return;
}



