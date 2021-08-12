
#include "TreeNode.h"
#include "MpapTime.h"
#include "TimeFunction.h"
#include "Functions.h"
#include "SolutionData.h"
#include "BasisFunctionsBSpline.h"
#include "myDataIntegrateCutFEM.h"
#include "myPoly.h"

#include "stabilisationRoutines.h"


extern  MpapTime  mpapTime;
extern List<TimeFunction> timeFunction;


/*
template<>
void TreeNode<3>::calcStiffnessAndResidualCutFEMFluid(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
    // GFEM for Stokes
    // with proper stabilisation
    //
    ///////////////////////////////////////

    int ii, jj, gp, nGauss, tempId;
    int TI, TIp1, TIp2, TIp3, TJ, TJp1, TJp2, TJp3;

    double  JacTemp, Jac, dvol, xx, yy, zz, stabParam;
    double  fact, fact2, b1, b2, b3, b4, b5, b6, b7, b8, acceFact;
    double  pres, Da, Db, af, am, muTaf, rad, urdr, urdr2, h2, tau[3];

    VectorXd  NN(totnlbf), dNN_dx(totnlbf), d2NN_dx2(totnlbf), dNN_dy(totnlbf), d2NN_dy2(totnlbf);
    VectorXd  dNN_dz(totnlbf), d2NN_dz2(totnlbf);
    VectorXd  N, dN_dx, d2N_dx2, dN_dy, d2N_dy2, dN_dz, d2N_dz2, d2N;
    VectorXd  res(4), res2(3), dp(3), Du(3), vel(3), velDot(3), force(3), gradTvel(3), rStab(4);
    MatrixXd  Dj(3, 4), grad(3,3), gradN(3,3), stress(3,3);
    myPoint  param, geom;
    Dj.setZero();

    double  rho = elmDat[3];
    double  mu  = elmDat[4];
    bforce[0]   = elmDat[5];
    bforce[1]   = elmDat[6];
    bforce[2]   = elmDat[7];

    af = SolnData->td(2);
    am = SolnData->td(1);
    acceFact = am*SolnData->td(9);
    dt = mpapTime.dt;

    double *gws;
    myPoint *gps;

    double  hx = bbox.maxBB[0]-bbox.minBB[0];
    double  hy = bbox.maxBB[1]-bbox.minBB[1];
    double  hz = bbox.maxBB[2]-bbox.minBB[2];

    volume = hx*hy*hz;

    matG.setZero();
    matG(0,0) = 4.0/hx/hx;
    matG(1,1) = 4.0/hy/hy;
    matG(2,2) = 4.0/hz/hz;

    if(domNums.size() > 1) // cut-cell
    {
      nGauss = Quadrature.gausspoints.size();
      
      gps = &(Quadrature.gausspoints[0]);
      gws = &(Quadrature.gaussweights[0]);
      
      //tempId = domainCur;
      JacTemp = 1.0;

      volume = 0.0;
      for(gp=0; gp<nGauss; gp++)
        volume += gws[gp];

      // For 2D problem
      // fact = sqrt(Vc/V); and fact = fact*fact;  ---->  fact = Vc/V;
      // For 3D problem
      // fact = (Vc/V)^(1/3). So, do fact = fact*fact;

      fact = volume/(hx*hy*hz);
      fact = pow(fact, 1.0/3.0);
      fact = fact*fact;

      matG(0,0) /= fact;
      matG(1,1) /= fact;
      matG(2,2) /= fact;
    }
    else
    {
      nGauss = GeomData->gausspoints.size();

      gps = &(GeomData->gausspoints[0]);
      gws = &(GeomData->gaussweights[0]);
      
      //tempId = domNums[0];
      JacTemp = JacMultElem;
    }

    //cout << " volume = " << volume << endl;
    //h2 = 3.0*volume/PI/4.0;
    //h2 = 4.0*volume/PI;

    //stabParam = h2/(12.0*mu)/degree[0]/degree[1]/degree[2];
    stabParam = h2/(12.0*mu);
    //cout << " stabParam = " << stabParam << endl;

    tau[0] = 0.0;
    tau[1] = elmDat[9]*stabParam;//rho;  // PSPG
    tau[2] = elmDat[10]*stabParam;//rho; // LSIC

    //cout << " nGauss " << nGauss << '\t' << tempId << endl;

    for(gp=0; gp<nGauss; gp++)
    {
        param[0]  = 0.5*(knotIncr[0] * gps[gp][0] + knotSum[0]);
        param[1]  = 0.5*(knotIncr[1] * gps[gp][1] + knotSum[1]);
        param[2]  = 0.5*(knotIncr[2] * gps[gp][2] + knotSum[2]);

        dvol = gws[gp] * JacTemp;

        geom[0] = GeomData->computeCoord(0, param[0]);
        geom[1] = GeomData->computeCoord(1, param[1]);
        geom[2] = GeomData->computeCoord(2, param[2]);

        //cout << uu << '\t' << vv << endl;
        //cout << xx << '\t' << yy << endl;

        //cout << gp << '\t' << gps[gp][0] << '\t' << gps[gp][1] << '\t' << gws[gp] << '\t' << dvol << endl;

        GeomData->computeBasisFunctions3D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy, dNN_dz, d2NN_dx2, d2NN_dy2, d2NN_dz2);

        if(parent == NULL)
        {
          //cout << " parent is NULL " << endl;
          N = NN;

          dN_dx = dNN_dx;
          dN_dy = dNN_dy;
          dN_dz = dNN_dz;

          d2N_dx2 = d2NN_dx2;
          d2N_dy2 = d2NN_dy2;
          d2N_dz2 = d2NN_dz2;
        }
        else
        {
          //cout << " parent is not NULL " << endl;
          N = SubDivMat*NN;

          dN_dx = SubDivMat*dNN_dx;
          dN_dy = SubDivMat*dNN_dy;
          dN_dz = SubDivMat*dNN_dz;

          d2N_dx2 = SubDivMat*d2NN_dx2;
          d2N_dy2 = SubDivMat*d2NN_dy2;
          d2N_dz2 = SubDivMat*d2NN_dz2;
        }

        //cout <<  " gpDomainId = " << gpDomainId << '\t' << totnlbf2 << endl;
        //printVector(GlobalBasisFuncs);

        d2N = d2N_dx2 + d2N_dy2 + d2N_dz2;

        vel.setZero();
        velDot.setZero();
        grad.setZero();
        dp.setZero();
        Du.setZero();
        pres = 0.0;

          for(ii=0;ii<totnlbf2;ii++)
          {
            TI   = ndof*GlobalBasisFuncs[ii];
            TIp1 = TI+1;
            TIp2 = TI+2;
            TIp3 = TI+3;

            b1     = SolnData->var1Cur(TI);
            b2     = SolnData->var1Cur(TIp1);
            b3     = SolnData->var1Cur(TIp2);
            b4     = SolnData->var1(TIp3);

            vel(0) += ( b1 * N[ii] );
            vel(1) += ( b2 * N[ii] );
            vel(2) += ( b3 * N[ii] );

            grad(0,0) += ( b1 * dN_dx[ii] );
            grad(0,1) += ( b1 * dN_dy[ii] );
            grad(0,2) += ( b1 * dN_dz[ii] );

            grad(1,0) += ( b2 * dN_dx[ii] );
            grad(1,1) += ( b2 * dN_dy[ii] );
            grad(1,2) += ( b2 * dN_dz[ii] );

            grad(2,0) += ( b3 * dN_dx[ii] );
            grad(2,1) += ( b3 * dN_dy[ii] );
            grad(2,2) += ( b3 * dN_dz[ii] );

            Du(0)  += ( b1 * d2N[ii] );
            Du(1)  += ( b2 * d2N[ii] );
            Du(2)  += ( b3 * d2N[ii] );

            pres   += ( b4 * N[ii] );

            dp(0)  += ( b4 * dN_dx[ii] );
            dp(1)  += ( b4 * dN_dy[ii] );
            dp(2)  += ( b4 * dN_dz[ii] );

            b1     = SolnData->var1DotCur(TI);
            b2     = SolnData->var1DotCur(TIp1);
            b3     = SolnData->var1DotCur(TIp2);

            velDot(0) += ( b1 * N[ii] );
            velDot(1) += ( b2 * N[ii] );
            velDot(2) += ( b3 * N[ii] );
          }

          // this is pseudo-stress
          stress = mu*grad;
          stress(0,0) -= pres;
          stress(1,1) -= pres;
          stress(2,2) -= pres;

          force.setZero();

          //force(0) = GeomData->analyDBC->computeForce(0, xx, yy);
          //force(1) = GeomData->analyDBC->computeForce(0, xx, yy);
          //force = 1.0;
          //force = 0.0;
          //cout << force(0) << '\t' << force(1) << endl;

          res2(0) = rho*( velDot(0) - force(0) ) ;
          res2(1) = rho*( velDot(1) - force(1) ) ;
          res2(2) = rho*( velDot(2) - force(2) ) ;

          rStab(0) = res2(0) - mu*Du(0) + dp(0) ;
          rStab(1) = res2(1) - mu*Du(1) + dp(1) ;
          rStab(2) = res2(2) - mu*Du(2) + dp(2) ;

          muTaf = mu*af;

          for(ii=0;ii<totnlbf2;ii++)
          {
            TI   = ndof*ii;
            TIp1 = TI+1;
            TIp2 = TI+2;
            TIp3 = TI+3;

            b1 = dN_dx[ii]*dvol;
            b2 = dN_dy[ii]*dvol;
            b3 = dN_dz[ii]*dvol;
            b4 = N[ii]*dvol;

            b5 = muTaf*b1;
            b6 = muTaf*b2;
            b7 = muTaf*b3;
            b8 = af*b4;

            for(jj=0;jj<totnlbf2;jj++)
            {
              TJ   = ndof*jj;
              TJp1 = TJ+1;
              TJp2 = TJ+2;
              TJp3 = TJ+3;

              fact2 = rho*acceFact*N[jj];

              // time acceleration term
              fact = b4*fact2 ;

              // diffusion term
              fact += (b5*dN_dx[jj] + b6*dN_dy[jj] + b7*dN_dz[jj]);

              Klocal(TI,   TJ)   += fact;
              Klocal(TIp1, TJp1) += fact;
              Klocal(TIp2, TJp2) += fact;

              gradN.setZero();

              // pressure term
              Klocal(TI,   TJp3) -= (b1*N[jj]);
              Klocal(TIp1, TJp3) -= (b2*N[jj]);
              Klocal(TIp2, TJp3) -= (b3*N[jj]);

              // continuity equation
              Klocal(TIp3, TJ)   += (b8*dN_dx[jj]);
              Klocal(TIp3, TJp1) += (b8*dN_dy[jj]);
              Klocal(TIp3, TJp2) += (b8*dN_dz[jj]);

              // PSPG stabilisation terms
              fact2 -= muTaf*d2N[jj];

              Dj(0,0) = gradN(0,0) + fact2;
              Dj(0,1) = gradN(0,1);
              Dj(0,2) = gradN(0,2);
              Dj(0,3) = dN_dx[jj];

              Dj(1,0) = gradN(1,0);
              Dj(1,1) = gradN(1,1) + fact2;
              Dj(1,2) = gradN(1,2);
              Dj(1,3) = dN_dy[jj];

              Dj(2,0) = gradN(2,0);
              Dj(2,1) = gradN(2,1);
              Dj(2,2) = gradN(2,2) + fact2;
              Dj(2,3) = dN_dz[jj];

              // PSPG
              Klocal(TIp3, TJ)   += ( b1*Dj(0,0) + b2*Dj(1,0) + b3*Dj(2,0) )*tau[1];
              Klocal(TIp3, TJp1) += ( b1*Dj(0,1) + b2*Dj(1,1) + b3*Dj(2,1) )*tau[1];
              Klocal(TIp3, TJp2) += ( b1*Dj(0,2) + b2*Dj(1,2) + b3*Dj(2,2) )*tau[1];
              Klocal(TIp3, TJp3) += ( b1*Dj(0,3) + b2*Dj(1,3) + b3*Dj(2,3) )*tau[1];

              // LSIC stabilisation

              Klocal(TI,   TJ)   += (b1*dN_dx[jj])*tau[2];
              Klocal(TI,   TJp1) += (b1*dN_dy[jj])*tau[2];
              Klocal(TI,   TJp2) += (b1*dN_dz[jj])*tau[2];

              Klocal(TIp1, TJ)   += (b2*dN_dx[jj])*tau[2];
              Klocal(TIp1, TJp1) += (b2*dN_dy[jj])*tau[2];
              Klocal(TIp1, TJp2) += (b2*dN_dz[jj])*tau[2];

              Klocal(TIp2, TJ)   += (b3*dN_dx[jj])*tau[2];
              Klocal(TIp2, TJp1) += (b3*dN_dy[jj])*tau[2];
              Klocal(TIp2, TJp2) += (b3*dN_dz[jj])*tau[2];
            }

            Flocal(TI)   -= (b4*res2(0) + b1*stress(0,0) + b2*stress(0,1) + b3*stress(0,2) );
            Flocal(TIp1) -= (b4*res2(1) + b1*stress(1,0) + b2*stress(1,1) + b3*stress(1,2) );
            Flocal(TIp2) -= (b4*res2(2) + b1*stress(2,0) + b2*stress(2,1) + b3*stress(2,2) );
            Flocal(TIp3) -= (b4*grad.trace());

            // PSPG stabilisation terms
            Flocal(TIp3) -= (tau[1]*( b1*rStab(0)+b2*rStab(1)+b3*rStab(2) ));

            // LSIC stabilisation terms
            Flocal(TI)   -= ( tau[2]*b1*grad.trace() );
            Flocal(TIp1) -= ( tau[2]*b2*grad.trace() );
            Flocal(TIp2) -= ( tau[2]*b3*grad.trace() );

          } // for(ii=0;ii<totnlbf2;ii++)
    }//gp

    return;
}
*/



/*
template<>
void TreeNode<3>::calcStiffnessAndResidualCutFEMFluid(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
    // GFEM for Navier-Stokes
    // with proper stabilisation
    // fully-implicit formulation
    ///////////////////////////////////////

    int  ii=0, jj=0, gp=0, nGauss=0;
    int  TI=0, TIp1=0, TIp2=0, TIp3=0, TJ=0, TJp1=0, TJp2=0, TJp3=0;

    double  JacTemp=0.0, Jac=0.0, dvol=0.0, c1=0.0;
    double  fact=0.0, fact2=0.0, b1=0.0, b2=0.0, b3=0.0, b4=0.0, b5=0.0, b6=0.0, b7=0.0, b8=0.0;
    double  pres=0.0, Da=0.0, Db=0.0, tau[3], beta[10];
    //double  CI=12.8435535945759;
    double  CI=10.0;

    VectorXd  NN(totnlbf), dNN_dx(totnlbf), d2NN_dx2(totnlbf), dNN_dy(totnlbf), d2NN_dy2(totnlbf);
    VectorXd  dNN_dz(totnlbf), d2NN_dz2(totnlbf);
    VectorXd  N, dN_dx, d2N_dx2, dN_dy, d2N_dy2, dN_dz, d2N_dz2, d2N, velTemp(3);
    VectorXd  res(4), res2(3), dp(3), Du(3), vel(3), velDot(3), force(3), gradTvel(3), rStab(4), velPrev(3);
    MatrixXd  Dj(3, 4), grad(3,3), gradN(3,3), stress(3,3);
    MatrixXd  matJ(3,3), matJinv(3,3), matG(3,3);
    myPoint  param, geom;
    Dj.setZero();

    double  rho = elmDat[3];
    double  mu  = elmDat[4];
    double  af = SolnData->td(2);
    double  am = SolnData->td(1);
    double  acceFact = am*SolnData->td(9);
    double  dt = mpapTime.dt;
    double  muTaf = mu*af;

    double *gws;
    myPoint *gps;

    double  hx = bbox.maxBB[0]-bbox.minBB[0];
    double  hy = bbox.maxBB[1]-bbox.minBB[1];
    double  hz = bbox.maxBB[2]-bbox.minBB[2];

    double  volume = hx*hy*hz;

    matG.setZero();
    matG(0,0) = 4.0/hx/hx;
    matG(1,1) = 4.0/hy/hy;
    matG(2,2) = 4.0/hz/hz;

    if(domNums.size() > 1) // cut-cell
    {
      nGauss = Quadrature.gausspoints.size();
      
      gps = &(Quadrature.gausspoints[0]);
      gws = &(Quadrature.gaussweights[0]);
      
      //tempId = domainCur;
      JacTemp = 1.0;

      volume = 0.0;
      for(gp=0; gp<nGauss; gp++)
        volume += gws[gp];

      // For 2D problem
      // fact = sqrt(Vc/V); and fact = fact*fact;  ---->  fact = Vc/V;
      // For 3D problem
      // fact = (Vc/V)^(1/3). So, do fact = fact*fact;

      fact = volume/(hx*hy*hz);
      fact = pow(fact, 1.0/3.0);
      fact = fact*fact;

      matG(0,0) /= fact;
      matG(1,1) /= fact;
      matG(2,2) /= fact;
    }
    else
    {
      nGauss = GeomData->gausspoints.size();

      gps = &(GeomData->gausspoints[0]);
      gws = &(GeomData->gaussweights[0]);
      
      //tempId = domNums[0];
      JacTemp = JacMultElem;
    }

    //cout << " nGauss " << nGauss << '\t' << tempId << endl;

    for(gp=0; gp<nGauss; gp++)
    {
        param[0]  = 0.5*(knotIncr[0] * gps[gp][0] + knotSum[0]);
        param[1]  = 0.5*(knotIncr[1] * gps[gp][1] + knotSum[1]);
        param[2]  = 0.5*(knotIncr[2] * gps[gp][2] + knotSum[2]);


        dvol = gws[gp] * JacTemp;

        geom[0] = GeomData->computeCoord(0, param[0]);
        geom[1] = GeomData->computeCoord(1, param[1]);
        geom[2] = GeomData->computeCoord(2, param[2]);

        //cout << uu << '\t' << vv << endl;
        //cout << xx << '\t' << yy << endl;

        //cout << gp << '\t' << gps[gp][0] << '\t' << gps[gp][1] << '\t' << gws[gp] << '\t' << dvol << endl;

        GeomData->computeBasisFunctions3D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy, dNN_dz, d2NN_dx2, d2NN_dy2, d2NN_dz2);

        if(parent == NULL)
        {
          //cout << " parent is NULL " << endl;
          N = NN;

          dN_dx = dNN_dx;
          dN_dy = dNN_dy;
          dN_dz = dNN_dz;

          d2N_dx2 = d2NN_dx2;
          d2N_dy2 = d2NN_dy2;
          d2N_dz2 = d2NN_dz2;
        }
        else
        {
          //cout << " parent is not NULL " << endl;
          N = SubDivMat*NN;

          dN_dx = SubDivMat*dNN_dx;
          dN_dy = SubDivMat*dNN_dy;
          dN_dz = SubDivMat*dNN_dz;

          d2N_dx2 = SubDivMat*d2NN_dx2;
          d2N_dy2 = SubDivMat*d2NN_dy2;
          d2N_dz2 = SubDivMat*d2NN_dz2;
        }

        //cout <<  " gpDomainId = " << gpDomainId << '\t' << totnlbf2 << endl;
        //printVector(GlobalBasisFuncs);

        d2N = d2N_dx2 + d2N_dy2 + d2N_dz2;

        vel.setZero();
        velPrev.setZero();
        velDot.setZero();
        grad.setZero();
        dp.setZero();
        Du.setZero();
        pres = 0.0;

          for(ii=0;ii<totnlbf2;ii++)
          {
            TI   = ndof*GlobalBasisFuncs[ii];
            TIp1 = TI+1;
            TIp2 = TI+2;
            TIp3 = TI+3;

            b1     = SolnData->var1Prev(TI);
            b2     = SolnData->var1Prev(TIp1);
            b3     = SolnData->var1Prev(TIp2);

            velPrev(0) += ( b1 * N[ii] );
            velPrev(1) += ( b2 * N[ii] );
            velPrev(2) += ( b3 * N[ii] );

            b1     = SolnData->var1Cur(TI);
            b2     = SolnData->var1Cur(TIp1);
            b3     = SolnData->var1Cur(TIp2);
            b4     = SolnData->var1(TIp3);

            vel(0) += ( b1 * N[ii] );
            vel(1) += ( b2 * N[ii] );
            vel(2) += ( b3 * N[ii] );

            grad(0,0) += ( b1 * dN_dx[ii] );
            grad(0,1) += ( b1 * dN_dy[ii] );
            grad(0,2) += ( b1 * dN_dz[ii] );

            grad(1,0) += ( b2 * dN_dx[ii] );
            grad(1,1) += ( b2 * dN_dy[ii] );
            grad(1,2) += ( b2 * dN_dz[ii] );

            grad(2,0) += ( b3 * dN_dx[ii] );
            grad(2,1) += ( b3 * dN_dy[ii] );
            grad(2,2) += ( b3 * dN_dz[ii] );

            Du(0)  += ( b1 * d2N[ii] );
            Du(1)  += ( b2 * d2N[ii] );
            Du(2)  += ( b3 * d2N[ii] );

            pres   += ( b4 * N[ii] );

            dp(0)  += ( b4 * dN_dx[ii] );
            dp(1)  += ( b4 * dN_dy[ii] );
            dp(2)  += ( b4 * dN_dz[ii] );

            b1     = SolnData->var1DotCur(TI);
            b2     = SolnData->var1DotCur(TIp1);
            b3     = SolnData->var1DotCur(TIp2);

            velDot(0) += ( b1 * N[ii] );
            velDot(1) += ( b2 * N[ii] );
            velDot(2) += ( b3 * N[ii] );
          }

          // this is pseudo-stress
          stress = mu*grad;
          stress(0,0) -= pres;
          stress(1,1) -= pres;
          stress(2,2) -= pres;

          force.setZero();

          //force(0) = GeomData->analyDBC->computeForce(0, xx, yy);
          //force(1) = GeomData->analyDBC->computeForce(0, xx, yy);
          //force = 1.0;
          //force = 0.0;
          //cout << force(0) << '\t' << force(1) << endl;

          gradTvel = grad*vel ;

          res2(0) = rho*( velDot(0) + gradTvel(0) - force(0) ) ;
          res2(1) = rho*( velDot(1) + gradTvel(1) - force(1) ) ;
          res2(2) = rho*( velDot(2) + gradTvel(2) - force(2) ) ;

          rStab(0) = res2(0) - mu*Du(0) + dp(0) ;
          rStab(1) = res2(1) - mu*Du(1) + dp(1) ;
          rStab(2) = res2(2) - mu*Du(2) + dp(2) ;

          // evaluate stabilisation parameters
          //
          velTemp(0) = velPrev(0);
          velTemp(1) = velPrev(1);
          velTemp(2) = velPrev(2);

          //velTemp(0) = vel(0);
          //velTemp(1) = vel(1);
          //velTemp(2) = vel(2);

          //evaluateStabParams_algo1(&velTemp(0), h, rho, mu, dt,  beta, tau);

          //evaluateStabParams_algo2(&velTemp(0), h, rho, mu, dt,  beta, tau);

          evaluateStabParams_algo3(velTemp, matG, dt, rho, mu, CI, tau);

          tau[0] *= elmDat[8];  // SUPG
          tau[1] *= elmDat[9];  // PSPG
          tau[2] *= elmDat[10]; // LSIC

          //cout << tau[0] << '\t' << tau[1] << '\t' << tau[2] << endl;

          for(ii=0;ii<totnlbf2;ii++)
          {
            TI   = ndof*ii;
            TIp1 = TI+1;
            TIp2 = TI+2;
            TIp3 = TI+3;

            b1 = dN_dx[ii]*dvol;
            b2 = dN_dy[ii]*dvol;
            b3 = dN_dz[ii]*dvol;
            b4 = N[ii]*dvol;

            b5 = muTaf*b1;
            b6 = muTaf*b2;
            b7 = muTaf*b3;
            b8 = af*b4;

            Da = (vel(0)*b1 + vel(1)*b2 + vel(2)*b3)*tau[0]*rho;

            for(jj=0;jj<totnlbf2;jj++)
            {
              TJ   = ndof*jj;
              TJp1 = TJ+1;
              TJp2 = TJ+2;
              TJp3 = TJ+3;

              fact2 = rho*acceFact*N[jj];

              // time acceleration term
              fact = b4*fact2 ;

              // diffusion term
              fact += (b5*dN_dx[jj] + b6*dN_dy[jj] + b7*dN_dz[jj]);

              Klocal(TI,   TJ)   += fact;
              Klocal(TIp1, TJp1) += fact;
              Klocal(TIp2, TJp2) += fact;

              // convection term

              gradN = grad*(rho*N[jj]);

              Db = rho*(vel(0)*dN_dx[jj] + vel(1)*dN_dy[jj] + vel(2)*dN_dz[jj] );

              gradN(0,0) += Db;
              gradN(1,1) += Db;
              gradN(2,2) += Db;
              gradN *= af;

              Klocal(TI,   TJ)   += (b4*gradN(0,0));
              Klocal(TI,   TJp1) += (b4*gradN(0,1));
              Klocal(TI,   TJp2) += (b4*gradN(0,2));

              Klocal(TIp1, TJ)   += (b4*gradN(1,0));
              Klocal(TIp1, TJp1) += (b4*gradN(1,1));
              Klocal(TIp1, TJp2) += (b4*gradN(1,2));

              Klocal(TIp2, TJ)   += (b4*gradN(2,0));
              Klocal(TIp2, TJp1) += (b4*gradN(2,1));
              Klocal(TIp2, TJp2) += (b4*gradN(2,2));

              // pressure term
              Klocal(TI,   TJp3) -= (b1*N[jj]);
              Klocal(TIp1, TJp3) -= (b2*N[jj]);
              Klocal(TIp2, TJp3) -= (b3*N[jj]);

              // continuity equation
              Klocal(TIp3, TJ)   += (b8*dN_dx[jj]);
              Klocal(TIp3, TJp1) += (b8*dN_dy[jj]);
              Klocal(TIp3, TJp2) += (b8*dN_dz[jj]);

              // SUPG and PSPG stabilisation terms
              fact2 -= muTaf*d2N[jj];

              Dj(0,0) = gradN(0,0) + fact2;
              Dj(0,1) = gradN(0,1);
              Dj(0,2) = gradN(0,2);
              Dj(0,3) = dN_dx[jj];

              Dj(1,0) = gradN(1,0);
              Dj(1,1) = gradN(1,1) + fact2;
              Dj(1,2) = gradN(1,2);
              Dj(1,3) = dN_dy[jj];

              Dj(2,0) = gradN(2,0);
              Dj(2,1) = gradN(2,1);
              Dj(2,2) = gradN(2,2) + fact2;
              Dj(2,3) = dN_dz[jj];

              // SUPG
              Klocal(TI, TJ)     += Da*Dj(0,0);
              Klocal(TI, TJp1)   += Da*Dj(0,1);
              Klocal(TI, TJp2)   += Da*Dj(0,2);
              Klocal(TI, TJp3)   += Da*Dj(0,3);

              Klocal(TIp1, TJ)   += Da*Dj(1,0);
              Klocal(TIp1, TJp1) += Da*Dj(1,1);
              Klocal(TIp1, TJp2) += Da*Dj(1,2);
              Klocal(TIp1, TJp3) += Da*Dj(1,3);

              Klocal(TIp2, TJ)   += Da*Dj(2,0);
              Klocal(TIp2, TJp1) += Da*Dj(2,1);
              Klocal(TIp2, TJp2) += Da*Dj(2,2);
              Klocal(TIp2, TJp3) += Da*Dj(2,3);

              c1 = tau[0] * af * rho * N[jj];

              Klocal(TI,   TJ)   += ( c1 * b1 * rStab(0) );  
              Klocal(TI,   TJp1) += ( c1 * b2 * rStab(0) );
              Klocal(TI,   TJp2) += ( c1 * b3 * rStab(0) );

              Klocal(TIp1, TJ)   += ( c1 * b1 * rStab(1) );
              Klocal(TIp1, TJp1) += ( c1 * b2 * rStab(1) );
              Klocal(TIp1, TJp2) += ( c1 * b3 * rStab(1) );

              Klocal(TIp2, TJ)   += ( c1 * b1 * rStab(2) );
              Klocal(TIp2, TJp1) += ( c1 * b2 * rStab(2) );
              Klocal(TIp2, TJp2) += ( c1 * b3 * rStab(2) );

              // PSPG
              Klocal(TIp3, TJ)   += ( b1*Dj(0,0) + b2*Dj(1,0) + b3*Dj(2,0) )*tau[1];
              Klocal(TIp3, TJp1) += ( b1*Dj(0,1) + b2*Dj(1,1) + b3*Dj(2,1) )*tau[1];
              Klocal(TIp3, TJp2) += ( b1*Dj(0,2) + b2*Dj(1,2) + b3*Dj(2,2) )*tau[1];
              Klocal(TIp3, TJp3) += ( b1*Dj(0,3) + b2*Dj(1,3) + b3*Dj(2,3) )*tau[1];

              // LSIC stabilisation

              fact2 = rho*af*tau[2];

              Klocal(TI,   TJ)   += (b1*dN_dx[jj])*fact2;
              Klocal(TI,   TJp1) += (b1*dN_dy[jj])*fact2;
              Klocal(TI,   TJp2) += (b1*dN_dz[jj])*fact2;

              Klocal(TIp1, TJ)   += (b2*dN_dx[jj])*fact2;
              Klocal(TIp1, TJp1) += (b2*dN_dy[jj])*fact2;
              Klocal(TIp1, TJp2) += (b2*dN_dz[jj])*fact2;

              Klocal(TIp2, TJ)   += (b3*dN_dx[jj])*fact2;
              Klocal(TIp2, TJp1) += (b3*dN_dy[jj])*fact2;
              Klocal(TIp2, TJp2) += (b3*dN_dz[jj])*fact2;
            }

            Flocal(TI)   -= (b4*res2(0) + b1*stress(0,0) + b2*stress(0,1) + b3*stress(0,2) );
            Flocal(TIp1) -= (b4*res2(1) + b1*stress(1,0) + b2*stress(1,1) + b3*stress(1,2) );
            Flocal(TIp2) -= (b4*res2(2) + b1*stress(2,0) + b2*stress(2,1) + b3*stress(2,2) );
            Flocal(TIp3) -= (b4*grad.trace());

            // SUPG stabilisation terms
            Flocal(TI)   -= Da*rStab(0);
            Flocal(TIp1) -= Da*rStab(1);
            Flocal(TIp2) -= Da*rStab(2);

            // PSPG stabilisation terms
            Flocal(TIp3) -= (tau[1]*( b1*rStab(0)+b2*rStab(1)+b3*rStab(2) ));

            // LSIC stabilisation terms
            fact2 = rho*grad.trace()*tau[2];

            Flocal(TI)   -= ( b1*fact2 );
            Flocal(TIp1) -= ( b2*fact2 );
            Flocal(TIp2) -= ( b3*fact2 );

          } // for(ii=0;ii<totnlbf2;ii++)
    }//gp

    return;
}
*/



//
template<>
void TreeNode<3>::calcStiffnessAndResidualCutFEMFluid(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
    // GFEM for Navier-Stokes
    // with proper stabilisation
    // semi-implicit formulation - type B
    ///////////////////////////////////////

    int  ii=0, jj=0, gp=0, nGauss=0;
    int  TI=0, TIp1=0, TIp2=0, TIp3=0, TJ=0, TJp1=0, TJp2=0, TJp3=0;

    double  JacTemp=0.0, Jac=0.0, dvol=0.0;
    double  fact=0.0, fact2=0.0, b1=0.0, b2=0.0, b3=0.0, b4=0.0, b5=0.0, b6=0.0, b7=0.0, b8=0.0;
    double  pres=0.0, Da=0.0, Db=0.0, tau[3];
    //double  CI=12.8435535945759;
    double  CI=10.0;

    VectorXd  NN(totnlbf), dNN_dx(totnlbf), d2NN_dx2(totnlbf), dNN_dy(totnlbf), d2NN_dy2(totnlbf);
    VectorXd  dNN_dz(totnlbf), d2NN_dz2(totnlbf);
    VectorXd  N(totnlbf), dN_dx(totnlbf), d2N_dx2(totnlbf), dN_dy(totnlbf), d2N_dy2(totnlbf), dN_dz(totnlbf), d2N_dz2(totnlbf), d2N(totnlbf), velTemp(3);
    VectorXd  res(4), res2(3), dp(3), Du(3), vel(3), velDot(3), force(3), gradTvel(3), rStab(4), velPrev(3);
    MatrixXd  Dj(3, 4), grad(3,3), gradN(3,3), stress(3,3), gradPrev(3,3);
    MatrixXd  matG(3,3);
    myPoint  param, geom;
    Dj.setZero();

    double  rho = elmDat[3];
    double  mu  = elmDat[4];
    double  af = SolnData->td(2);
    double  am = SolnData->td(1);
    double  acceFact = am*SolnData->td(9);
    double  dt = mpapTime.dt;
    double  muTaf = mu*af;

    double *gws;
    myPoint *gps;

    double  hx = bbox.maxBB[0]-bbox.minBB[0];
    double  hy = bbox.maxBB[1]-bbox.minBB[1];
    double  hz = bbox.maxBB[2]-bbox.minBB[2];

    double  volume = hx*hy*hz;

    matG.setZero();
    matG(0,0) = 4.0/hx/hx;
    matG(1,1) = 4.0/hy/hy;
    matG(2,2) = 4.0/hz/hz;

    if(domNums.size() > 1) // cut-cell
    {
      nGauss = Quadrature.gausspoints.size();
      
      gps = &(Quadrature.gausspoints[0]);
      gws = &(Quadrature.gaussweights[0]);
      
      //tempId = domainCur;
      JacTemp = 1.0;

      volume = 0.0;
      for(gp=0; gp<nGauss; gp++)
        volume += gws[gp];

      // For 2D problem
      // fact = sqrt(Vc/V); and fact = fact*fact;  ---->  fact = Vc/V;
      // For 3D problem
      // fact = (Vc/V)^(1/3). So, do fact = fact*fact;

      fact = volume/(hx*hy*hz);
      fact = pow(fact, 1.0/3.0);
      fact = fact*fact;

      matG(0,0) /= fact;
      matG(1,1) /= fact;
      matG(2,2) /= fact;
    }
    else
    {
      nGauss = GeomData->gausspoints.size();

      gps = &(GeomData->gausspoints[0]);
      gws = &(GeomData->gaussweights[0]);

      //tempId = domNums[0];
      JacTemp = JacMultElem;
    }

    for(gp=0; gp<nGauss; gp++)
    {
        param[0]  = 0.5*(knotIncr[0] * gps[gp][0] + knotSum[0]);
        param[1]  = 0.5*(knotIncr[1] * gps[gp][1] + knotSum[1]);
        param[2]  = 0.5*(knotIncr[2] * gps[gp][2] + knotSum[2]);

        dvol = gws[gp] * JacTemp;

        geom[0] = GeomData->computeCoord(0, param[0]);
        geom[1] = GeomData->computeCoord(1, param[1]);
        geom[2] = GeomData->computeCoord(2, param[2]);

        if(parent == NULL) //cout << " parent is NULL " << endl;
        {
          GeomData->computeBasisFunctions3D(knotBegin, knotIncr, param, N, dN_dx, dN_dy, dN_dz, d2N_dx2, d2N_dy2, d2N_dz2);
        }
        else //cout << " parent is not NULL " << endl;
        {
          GeomData->computeBasisFunctions3D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy, dNN_dz, d2NN_dx2, d2NN_dy2, d2NN_dz2);

          N = SubDivMat*NN;

          dN_dx = SubDivMat*dNN_dx;
          dN_dy = SubDivMat*dNN_dy;
          dN_dz = SubDivMat*dNN_dz;

          d2N_dx2 = SubDivMat*d2NN_dx2;
          d2N_dy2 = SubDivMat*d2NN_dy2;
          d2N_dz2 = SubDivMat*d2NN_dz2;
        }

        d2N = d2N_dx2 + d2N_dy2 + d2N_dz2;

        vel.setZero();
        velPrev.setZero();
        velDot.setZero();
        grad.setZero();
        gradPrev.setZero();
        dp.setZero();
        Du.setZero();
        pres = 0.0;

          for(ii=0;ii<totnlbf2;ii++)
          {
            TI   = ndof*GlobalBasisFuncs[ii];
            TIp1 = TI+1;
            TIp2 = TI+2;
            TIp3 = TI+3;

            b1   = SolnData->var1Prev(TI);
            b2   = SolnData->var1Prev(TIp1);
            b3   = SolnData->var1Prev(TIp2);

            velPrev(0) += ( b1 * N[ii] );
            velPrev(1) += ( b2 * N[ii] );
            velPrev(2) += ( b3 * N[ii] );

            gradPrev(0,0) += ( b1 * dN_dx[ii] );
            gradPrev(0,1) += ( b1 * dN_dy[ii] );
            gradPrev(0,2) += ( b1 * dN_dz[ii] );

            gradPrev(1,0) += ( b2 * dN_dx[ii] );
            gradPrev(1,1) += ( b2 * dN_dy[ii] );
            gradPrev(1,2) += ( b2 * dN_dz[ii] );

            gradPrev(2,0) += ( b3 * dN_dx[ii] );
            gradPrev(2,1) += ( b3 * dN_dy[ii] );
            gradPrev(2,2) += ( b3 * dN_dz[ii] );

            b1     = SolnData->var1Cur(TI);
            b2     = SolnData->var1Cur(TIp1);
            b3     = SolnData->var1Cur(TIp2);
            b4     = SolnData->var1(TIp3);

            vel(0) += ( b1 * N[ii] );
            vel(1) += ( b2 * N[ii] );
            vel(2) += ( b3 * N[ii] );

            grad(0,0) += ( b1 * dN_dx[ii] );
            grad(0,1) += ( b1 * dN_dy[ii] );
            grad(0,2) += ( b1 * dN_dz[ii] );

            grad(1,0) += ( b2 * dN_dx[ii] );
            grad(1,1) += ( b2 * dN_dy[ii] );
            grad(1,2) += ( b2 * dN_dz[ii] );

            grad(2,0) += ( b3 * dN_dx[ii] );
            grad(2,1) += ( b3 * dN_dy[ii] );
            grad(2,2) += ( b3 * dN_dz[ii] );

            Du(0)  += ( b1 * d2N[ii] );
            Du(1)  += ( b2 * d2N[ii] );
            Du(2)  += ( b3 * d2N[ii] );

            pres   += ( b4 * N[ii] );

            dp(0)  += ( b4 * dN_dx[ii] );
            dp(1)  += ( b4 * dN_dy[ii] );
            dp(2)  += ( b4 * dN_dz[ii] );

            b1     = SolnData->var1DotCur(TI);
            b2     = SolnData->var1DotCur(TIp1);
            b3     = SolnData->var1DotCur(TIp2);

            velDot(0) += ( b1 * N[ii] );
            velDot(1) += ( b2 * N[ii] );
            velDot(2) += ( b3 * N[ii] );
          }

          // this is pseudo-stress
          stress = mu*grad;
          stress(0,0) -= pres;
          stress(1,1) -= pres;
          stress(2,2) -= pres;

          force.setZero();

          //force(0) = GeomData->analyDBC->computeForce(0, xx, yy);
          //force(1) = GeomData->analyDBC->computeForce(0, xx, yy);
          //force = 1.0;
          //force = 0.0;
          //cout << force(0) << '\t' << force(1) << endl;

          gradTvel = gradPrev*vel + grad*velPrev - gradPrev*velPrev;

          res2(0) = rho*( velDot(0) + gradTvel(0) - force(0) ) ;
          res2(1) = rho*( velDot(1) + gradTvel(1) - force(1) ) ;
          res2(2) = rho*( velDot(2) + gradTvel(2) - force(2) ) ;

          rStab(0) = res2(0) - mu*Du(0) + dp(0) ;
          rStab(1) = res2(1) - mu*Du(1) + dp(1) ;
          rStab(2) = res2(2) - mu*Du(2) + dp(2) ;

          // evaluate stabilisation parameters
          //
          velTemp(0) = velPrev(0);
          velTemp(1) = velPrev(1);
          velTemp(2) = velPrev(2);

          //evaluateStabParams_algo1(&velTemp(0), h, rho, mu, dt,  beta, tau);
          //evaluateStabParams_algo2(&velTemp(0), h, rho, mu, dt,  beta, tau);
          evaluateStabParams_algo3(velTemp, matG, dt, rho, mu, CI, tau);

          tau[0] *= elmDat[8];  // SUPG
          tau[1] *= elmDat[9];  // PSPG
          tau[2] *= elmDat[10]; // LSIC

          //cout << tau[0] << '\t' << tau[1] << '\t' << tau[2] << endl;

          for(ii=0;ii<totnlbf2;ii++)
          {
            TI   = ndof*ii;
            TIp1 = TI+1;
            TIp2 = TI+2;
            TIp3 = TI+3;

            b1 = dN_dx[ii]*dvol;
            b2 = dN_dy[ii]*dvol;
            b3 = dN_dz[ii]*dvol;
            b4 = N[ii]*dvol;

            b5 = muTaf*b1;
            b6 = muTaf*b2;
            b7 = muTaf*b3;
            b8 = af*b4;

            Da = (velPrev(0)*b1 + velPrev(1)*b2 + velPrev(2)*b3)*tau[0]*rho;

            for(jj=0;jj<totnlbf2;jj++)
            {
              TJ   = ndof*jj;
              TJp1 = TJ+1;
              TJp2 = TJ+2;
              TJp3 = TJ+3;

              fact2 = rho*acceFact*N[jj];

              // time acceleration term
              fact = b4*fact2 ;

              // diffusion term
              fact += (b5*dN_dx[jj] + b6*dN_dy[jj] + b7*dN_dz[jj]);

              Klocal(TI,   TJ)   += fact;
              Klocal(TIp1, TJp1) += fact;
              Klocal(TIp2, TJp2) += fact;

              // convection term

              gradN = gradPrev*(rho*N[jj]);

              Db = rho*(velPrev(0)*dN_dx[jj] + velPrev(1)*dN_dy[jj] + velPrev(2)*dN_dz[jj] );

              gradN(0,0) += Db;
              gradN(1,1) += Db;
              gradN(2,2) += Db;

              Klocal(TI,   TJ)   += (b8*gradN(0,0));
              Klocal(TI,   TJp1) += (b8*gradN(0,1));
              Klocal(TI,   TJp2) += (b8*gradN(0,2));

              Klocal(TIp1, TJ)   += (b8*gradN(1,0));
              Klocal(TIp1, TJp1) += (b8*gradN(1,1));
              Klocal(TIp1, TJp2) += (b8*gradN(1,2));

              Klocal(TIp2, TJ)   += (b8*gradN(2,0));
              Klocal(TIp2, TJp1) += (b8*gradN(2,1));
              Klocal(TIp2, TJp2) += (b8*gradN(2,2));

              // pressure term
              Klocal(TI,   TJp3) -= (b1*N[jj]);
              Klocal(TIp1, TJp3) -= (b2*N[jj]);
              Klocal(TIp2, TJp3) -= (b3*N[jj]);

              // continuity equation
              Klocal(TIp3, TJ)   += (b8*dN_dx[jj]);
              Klocal(TIp3, TJp1) += (b8*dN_dy[jj]);
              Klocal(TIp3, TJp2) += (b8*dN_dz[jj]);

              // SUPG and PSPG stabilisation terms
              fact2 -= muTaf*d2N[jj];

              gradN *= af;

              Dj(0,0) = gradN(0,0) + fact2;
              Dj(0,1) = gradN(0,1);
              Dj(0,2) = gradN(0,2);
              Dj(0,3) = dN_dx[jj];

              Dj(1,0) = gradN(1,0);
              Dj(1,1) = gradN(1,1) + fact2;
              Dj(1,2) = gradN(1,2);
              Dj(1,3) = dN_dy[jj];

              Dj(2,0) = gradN(2,0);
              Dj(2,1) = gradN(2,1);
              Dj(2,2) = gradN(2,2) + fact2;
              Dj(2,3) = dN_dz[jj];

              // SUPG
              Klocal(TI, TJ)     += Da*Dj(0,0);
              Klocal(TI, TJp1)   += Da*Dj(0,1);
              Klocal(TI, TJp2)   += Da*Dj(0,2);
              Klocal(TI, TJp3)   += Da*Dj(0,3);

              Klocal(TIp1, TJ)   += Da*Dj(1,0);
              Klocal(TIp1, TJp1) += Da*Dj(1,1);
              Klocal(TIp1, TJp2) += Da*Dj(1,2);
              Klocal(TIp1, TJp3) += Da*Dj(1,3);

              Klocal(TIp2, TJ)   += Da*Dj(2,0);
              Klocal(TIp2, TJp1) += Da*Dj(2,1);
              Klocal(TIp2, TJp2) += Da*Dj(2,2);
              Klocal(TIp2, TJp3) += Da*Dj(2,3);

              // PSPG
              Klocal(TIp3, TJ)   += ( b1*Dj(0,0) + b2*Dj(1,0) + b3*Dj(2,0) )*tau[1];
              Klocal(TIp3, TJp1) += ( b1*Dj(0,1) + b2*Dj(1,1) + b3*Dj(2,1) )*tau[1];
              Klocal(TIp3, TJp2) += ( b1*Dj(0,2) + b2*Dj(1,2) + b3*Dj(2,2) )*tau[1];
              Klocal(TIp3, TJp3) += ( b1*Dj(0,3) + b2*Dj(1,3) + b3*Dj(2,3) )*tau[1];

              // LSIC stabilisation

              fact2 = rho*af*tau[2];

              Klocal(TI,   TJ)   += (b1*dN_dx[jj])*fact2;
              Klocal(TI,   TJp1) += (b1*dN_dy[jj])*fact2;
              Klocal(TI,   TJp2) += (b1*dN_dz[jj])*fact2;

              Klocal(TIp1, TJ)   += (b2*dN_dx[jj])*fact2;
              Klocal(TIp1, TJp1) += (b2*dN_dy[jj])*fact2;
              Klocal(TIp1, TJp2) += (b2*dN_dz[jj])*fact2;

              Klocal(TIp2, TJ)   += (b3*dN_dx[jj])*fact2;
              Klocal(TIp2, TJp1) += (b3*dN_dy[jj])*fact2;
              Klocal(TIp2, TJp2) += (b3*dN_dz[jj])*fact2;
            }

            Flocal(TI)   -= (b4*res2(0) + b1*stress(0,0) + b2*stress(0,1) + b3*stress(0,2) );
            Flocal(TIp1) -= (b4*res2(1) + b1*stress(1,0) + b2*stress(1,1) + b3*stress(1,2) );
            Flocal(TIp2) -= (b4*res2(2) + b1*stress(2,0) + b2*stress(2,1) + b3*stress(2,2) );
            Flocal(TIp3) -= (b4*grad.trace());

            // SUPG stabilisation terms
            Flocal(TI)   -= Da*rStab(0);
            Flocal(TIp1) -= Da*rStab(1);
            Flocal(TIp2) -= Da*rStab(2);

            // PSPG stabilisation terms
            Flocal(TIp3) -= (tau[1]*( b1*rStab(0)+b2*rStab(1)+b3*rStab(2) ));

            // LSIC stabilisation terms
            fact2 = rho*grad.trace()*tau[2];

            Flocal(TI)   -= ( b1*fact2 );
            Flocal(TIp1) -= ( b2*fact2 );
            Flocal(TIp2) -= ( b3*fact2 );

          } // for(ii=0;ii<totnlbf2;ii++)
    }//gp

    return;
}
//




template<>
void TreeNode<3>::applyBoundaryConditionsAtApointCutFEMFluid(myDataIntegrateCutFEM& myData)
{
  // compute stiffness and force vectors corresponding 
  // to Nitsche method of applying interface conditions
  // 
  // diagonal terms
  
    int ii, jj, TI, TIp1, TIp2, TIp3, TJ, TJp1, TJp2, TJp3;

    double  fact, fact1, fact2, fact3, res, trac;
    double  Ta1, Ta2, Ta3, Tb1, Tb2, Tb3, mu1, mu2, specVal;
    double  bb1, bb2, u1, u2, t1, t2, pres;

    VectorXd  NN(totnlbf), dNN_dx(totnlbf), dNN_dy(totnlbf), dNN_dz(totnlbf);
    VectorXd  N, dN_dx, dN_dy, dN_dz;
    VectorXd  vel1(3), vel2(3), trac1(3), trac2(3);
    MatrixXd  stress1(3,3), stress2(3,3);
    myPoint   normal1, normal2;

    double  rho = elmDat[3];
    double  af = SolnData->td(2);

    GeomData->computeBasisFunctions2D(knotBegin, knotIncr, myData.param, NN, dNN_dx, dNN_dy, dNN_dz);

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
    
    // at this moment no jumps in velocity or tractions is considered
    
    vel1[0] = myData.specVal[0];
    vel1[1] = myData.specVal[1];
    vel1[2] = myData.specVal[2];
    vel2 = vel1;

    vel1(0) -= computeValueCur(0, N);
    vel1(1) -= computeValueCur(1, N);
    vel1(2) -= computeValueCur(2, N);

    vel2(0) -= computeValue2Cur(0, N);
    vel2(1) -= computeValue2Cur(1, N);
    vel2(2) -= computeValue2Cur(2, N);

    ///////////////////////////////////
    // compute the stresses and tractions
    ///////////////////////////////////

    // domain 1

    mu1 = elmDat[4];
    normal1 = -myData.normal;

    stress1(0,0) = mu1*computeValueCur(0, dN_dx);
    stress1(0,1) = mu1*computeValueCur(0, dN_dy);
    stress1(0,2) = mu1*computeValueCur(0, dN_dz);

    stress1(1,0) = mu1*computeValueCur(1, dN_dx);
    stress1(1,1) = mu1*computeValueCur(1, dN_dy);
    stress1(1,2) = mu1*computeValueCur(1, dN_dz);

    stress1(2,0) = mu1*computeValueCur(2, dN_dx);
    stress1(2,1) = mu1*computeValueCur(2, dN_dy);
    stress1(2,2) = mu1*computeValueCur(2, dN_dz);

    pres         = computeValue(2, N);

    stress1(0,0) -= pres;
    stress1(1,1) -= pres;
    stress1(2,2) -= pres;

    trac1[0] = 0.0 + ( stress1(0,0)*normal1[0] + stress1(0,1)*normal1[1] + stress1(0,2)*normal1[2] );
    trac1[1] = 0.0 + ( stress1(1,0)*normal1[0] + stress1(1,1)*normal1[1] + stress1(1,2)*normal1[2] );
    trac1[2] = 0.0 + ( stress1(2,0)*normal1[0] + stress1(2,1)*normal1[1] + stress1(2,2)*normal1[2] );

    // domain 2

    mu2 = elmDat[5];
    normal2 = -normal1;

    stress2(0,0) = mu2*computeValue2Cur(0, dN_dx);
    stress2(0,1) = mu2*computeValue2Cur(0, dN_dy);
    stress2(0,2) = mu2*computeValue2Cur(0, dN_dz);

    stress2(1,0) = mu2*computeValue2Cur(1, dN_dx);
    stress2(1,1) = mu2*computeValue2Cur(1, dN_dy);
    stress2(1,2) = mu2*computeValue2Cur(1, dN_dz);

    stress2(2,0) = mu2*computeValue2Cur(2, dN_dx);
    stress2(2,1) = mu2*computeValue2Cur(2, dN_dy);
    stress2(2,2) = mu2*computeValue2Cur(2, dN_dz);

    pres         = computeValue2(2, N);

    stress2(0,0) -= pres;
    stress2(1,1) -= pres;
    stress2(2,2) -= pres;

    trac2[0] = 0.0 + ( stress2(0,0)*normal2[0] + stress1(0,1)*normal2[1] + stress2(0,2)*normal2[2] );
    trac2[1] = 0.0 + ( stress2(1,0)*normal2[0] + stress1(1,1)*normal2[1] + stress2(1,2)*normal2[2] );
    trac2[2] = 0.0 + ( stress2(2,0)*normal2[0] + stress1(2,1)*normal2[1] + stress2(2,2)*normal2[2] );

    //for(ii=0;ii<totnlbf;ii++)
      //printf(" \t %14.8f \n", N[ii]);

    //specVal = GeomData->analyDBC->computeValue(0, myData.geom[0], myData.geom[1]);

    //cout << " myData.PENALTY = " << myData.PENALTY << endl;

    for(ii=0;ii<totnlbf2;ii++)
    {
        TI   = ndof*ii;
        TIp1 = TI+1;
        TIp2 = TI+2;
        TIp3 = TI+3;

        bb1 = N[ii] * myData.dvol;
        bb2 = bb1 * myData.PENALTY ;

        Ta1 = (myData.dvol*mu1)*( dN_dx[ii]*normal1[0] + dN_dy[ii]*normal1[1] + dN_dz[ii]*normal1[2]);

        for(jj=0;jj<totnlbf2;jj++)
        {
          TJ   = ndof*jj;
          TJp1 = TJ+1;
          TJp2 = TJ+2;
          TJp3 = TJ+3;

          fact = af* bb2 * N[jj];
          // stabilisation term
          myData.K1(TI,   TJ)   += fact;
          myData.K1(TIp1, TJp1) += fact;
          myData.K1(TIp2, TJp2) += fact;

          // Nitsche terms
          Tb1 = af*mu1*( dN_dx[jj]*normal1[0] + dN_dy[jj]*normal1[1] + dN_dz[jj]*normal1[2] );

          fact1 = -normal1[0]*N[jj];
          fact2 = -normal1[1]*N[jj];
          fact3 = -normal1[2]*N[jj];

          myData.K1(TI, TJ)      -= (bb1*Tb1);
          myData.K1(TI, TJp3)    -= (bb1*fact1);

          myData.K1(TI, TJ)      -= (Ta1*af*N[jj])*myData.NitscheFact;
          myData.K1(TIp3, TJ)    -= (bb1*fact1*af*myData.NitscheFact);

          myData.K1(TIp1, TJp1)  -= (bb1*Tb1);
          myData.K1(TIp1, TJp3)  -= (bb1*fact2);

          myData.K1(TIp1, TJp1)  -= (Ta1*af*N[jj])*myData.NitscheFact;
          myData.K1(TIp3, TJp1)  -= (bb1*fact2*af*myData.NitscheFact);

          myData.K1(TIp2, TJp2)  -= (bb1*Tb1);
          myData.K1(TIp2, TJp3)  -= (bb1*fact3);

          myData.K1(TIp2, TJp2)  -= (Ta1*af*N[jj])*myData.NitscheFact;
          myData.K1(TIp3, TJp2)  -= (bb1*fact3*af*myData.NitscheFact);
        }

        // stabilisation terms
        myData.F1(TI)   += (bb2*vel1[0]);
        myData.F1(TIp1) += (bb2*vel1[1]);
        myData.F1(TIp2) += (bb2*vel1[2]);

        // Nitsche terms

        myData.F1(TI)   -= (-bb1*trac1[0]);
        myData.F1(TI)   -= (Ta1*vel1[0])*myData.NitscheFact;
        myData.F1(TIp3) -= (bb1*(-normal1[0]*vel1[0]))*myData.NitscheFact;

        myData.F1(TIp1) -= (-bb1*trac1[1]);
        myData.F1(TIp1) -= (Ta1*vel1[1])*myData.NitscheFact;
        myData.F1(TIp3) -= (bb1*(-normal1[1]*vel1[1]))*myData.NitscheFact;

        myData.F1(TIp2) -= (-bb1*trac1[2]);
        myData.F1(TIp2) -= (Ta1*vel1[2])*myData.NitscheFact;
        myData.F1(TIp3) -= (bb1*(-normal1[2]*vel1[2]))*myData.NitscheFact;

        // terms for the monolithic scheme

        //myData.K2(0, TIp1)  -= (af*bb2);
        //myData.K2(0, TIp1)  += (af*Ta1);
        //myData.K2(0, TIp2)  -= (bb1*normal1[1]);

        //myData.Kc(0, 0)     += (af*myData.dvol*myData.PENALTY);

        //myData.F2(0)   -= (myData.dvol*myData.PENALTY*vel1[1]);
        //myData.F2(0)   -= (myData.dvol*trac1[1]);
    } // for(ii=0;ii<totnlbf2;ii++)
    //

   return;
}



template<>
void TreeNode<3>::applyGhostPenaltyCutFEM(myDataIntegrateCutFEM& myData)
{

   return;
}



template<>
void TreeNode<3>::applyBoundaryConditionsAtApointCutFEMFluid2(myDataIntegrateCutFEM& myData)
{
  // compute stiffness and force vectors corresponding 
  // to Nitsche method of applying interface conditions
  // 
  // coupling terms


   return;
}



template<>
int TreeNode<3>::computeGaussPointsAdapIntegrationBoundary(int side, int refLev1, int refLev2, int inclFlag, int domainCur)
{
    myPoint  pt1, pt2, ptTemp;

    int  domTemp, ii, jj, kk, gp, e1, e2, nGauss;

    double  val1, val2;
    bool  flag;
    
    AdaptiveBinarytree<2>  *AdapIntgCellBoundary = new AdaptiveBinarytree<2>(0);

    switch(side)
    {
        case 0:
        
          e1 = 1;  e2 = 2;
          
          val1 = knotBegin[0];
          val2 = GeomData->computeCoord(0, 0.0);

        break;

        case 1:

          e1 = 1;  e2 = 2;
          
          val1 = knotEnd[0];
          val2 = GeomData->computeCoord(0, 1.0);

        break;

        case 2:

          e1 = 0;  e2 = 2;

          val1 = knotBegin[1];
          val2 = GeomData->computeCoord(1, 0.0);

        break;

        case 3:

          e1 = 0;  e2 = 2;

          val1 = knotEnd[1];
          val2 = GeomData->computeCoord(1, 1.0);

        break;

        case 4:

          e1 = 0;  e2 = 1;

          val1 = knotBegin[2];
          val2 = GeomData->computeCoord(2, 0.0);

        break;

        case 5:

          e1 = 0;  e2 = 1;

          val1 = knotEnd[2];
          val2 = GeomData->computeCoord(2, 1.0);

        break;

        default :

          cout << " Invalid 'side' value in GeomDataHBSplines::getBoundaryGPs2D " << endl;
        break;
    } //switch(side)

    // identify the edge of the subtriangle that is on the boundary of the background grid

    AdapIntgCellBoundary->setKnots(knotBegin[e1], knotEnd[e1], knotBegin[e2], knotEnd[e2]);
    AdapIntgCellBoundary->setSideTemp(side);
    AdapIntgCellBoundary->setParam3(val1);
    AdapIntgCellBoundary->setCoord3(val2);

    AdapIntgCellBoundary->GeomData = GeomData;
    AdapIntgCellBoundary->domNums = domNums;
    AdapIntgCellBoundary->setSplitDirection(-1);

    //cout << " AAAAAAAAAA .... " << refLev1 << endl;
  
    AdapIntgCellBoundary->prepareData();
    AdapIntgCellBoundary->subDivide(refLev1);
  
    BoundaryQuadrature[side].reset();

    AdapIntgCellBoundary->computeGaussPoints2Dfor3D(refLev2, 0, 1, 1, BoundaryQuadrature[side]);

    // parametric domain to integration master-quadrilateral domain

    nGauss = BoundaryQuadrature[side].gausspoints.size();

    for(gp=0; gp<nGauss; gp++)
    {
      for(ii=0; ii<3; ii++)
        BoundaryQuadrature[side].gausspoints[gp][ii] = (2.0*BoundaryQuadrature[side].gausspoints[gp][ii] - knotSum[ii])/knotIncr[ii];
    }

    //printVector(boundaryGPs1);
    //printVector(boundaryGPs2);
//
    return 1;
}




template<>
void TreeNode<3>::applyDirichletBCsCutFEMFluid(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
  if(DirichletData.size() > 0)
  {
    int ii, jj, aa, gp, nGauss, TI, TIp1, TIp2, TIp3, index, TJ, TJp1, TJp2, TJp3, dir, side, profile_type;
    double  theta, ll1, ul1, ll2, ul2, Ta, Tb, res, JacTemp, NitscheFact, pres;
    double  dvol, specVal, PENALTY, Jac, fact, bb1, bb2, bb3;
    double  xc, yc, zc, R;
    bool  isNitsche;

    double  rho = elmDat[3];
    double  mu  = elmDat[4];
    //double  af = SolnData->td(2);
    double  af = 1.0;

    VectorXd  N(totnlbf), dN_dx(totnlbf), dN_dy(totnlbf), dN_dz(totnlbf);
    VectorXd  NN(totnlbf), dNN_dx(totnlbf), dNN_dy(totnlbf), dNN_dz(totnlbf), trac(3);
    MatrixXd  grad(3,3), stress(3,3);
    myPoint  param, geom, normal;
    double *gws;
    myPoint *gps;

    for(aa=0;aa<DirichletData.size();aa++)
    {
        // printVector(DirichletData[aa]);

        isNitsche   = false;
        side        = (int) (DirichletData[aa][0] - 1);
        dir         = (int) (DirichletData[aa][1] - 1);

        PENALTY     = DirichletData[aa][2];
        isNitsche   = ( (int) DirichletData[aa][3] == 1 );
        NitscheFact = DirichletData[aa][4];

        //PENALTY    = 1.0/max(hx, hy);                     // GeomData-> degree;

        //for symmetric Nitsche method   -> NitscheFact =  1.0
        //for unsymmetric Nitsche method -> NitscheFact = -1.0

        profile_type  = (int) DirichletData[aa][5];
        specVal       = DirichletData[aa][6];
        ll1           = DirichletData[aa][7];
        ul1           = DirichletData[aa][8];
        ll2           = DirichletData[aa][9];
        ul2           = DirichletData[aa][10];


        normal = GeomData->boundaryNormals[side];

        if(domNums.size() > 1)
        {
          nGauss = BoundaryQuadrature[side].gausspoints.size() ;

          gps = &(BoundaryQuadrature[side].gausspoints[0]);
          gws = &(BoundaryQuadrature[side].gaussweights[0]);

          JacTemp = 1.0;
        }
        else
        {
          nGauss = GeomData->boundaryQuadrature3D[side].gausspoints.size();

          gps = &(GeomData->boundaryQuadrature3D[side].gausspoints[0]);
          gws = &(GeomData->boundaryQuadrature3D[side].gaussweights[0]);

          JacTemp = GeomData->boundaryJacobians[side][level];
        }

        for(gp=0; gp<nGauss; gp++)
        {
            param[0] = 0.5*(knotIncr[0] * gps[gp][0] + knotSum[0]);
            param[1] = 0.5*(knotIncr[1] * gps[gp][1] + knotSum[1]);
            param[2] = 0.5*(knotIncr[2] * gps[gp][2] + knotSum[2]);

            dvol = JacTemp * gws[gp] ;

            geom[0] = GeomData->computeCoord(0, param[0]);
            geom[1] = GeomData->computeCoord(1, param[1]);
            geom[2] = GeomData->computeCoord(2, param[2]);

            //printf(" %4d \t %4d \t %12.6f \t %12.6f \t %12.6f \n", side, dir, geom[0], geom[1], geom[2]);

            //rad = sqrt((yy-yc)*(yy-yc)+(zz-zc)*(zz-zc));

            //rad = sqrt(xx*xx+yy*yy);
            //cout << xx << '\t' << yy << endl;

            //if(! GeomData->polyImm.within(xx, yy) )
            //{

            if(parent == NULL)
            {
              GeomData->computeBasisFunctions3D(knotBegin, knotIncr, param, N, dN_dx, dN_dy, dN_dz);
            }
            else
            {
              GeomData->computeBasisFunctions3D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy, dNN_dz);

              N = SubDivMat*NN;
              dN_dx = SubDivMat*dNN_dx;
              dN_dy = SubDivMat*dNN_dy;
              dN_dz = SubDivMat*dNN_dz;
            }

            // constant value
            if(profile_type == 1 )
            {
                specVal = DirichletData[aa][6];
            }
            // parabolic profile
            else
            {
                // faces perpendicular to X-axis
                if(side == 0 || side == 1)
                  specVal = DirichletData[aa][6]*(6.0/(ul1-ll1)/(ul1-ll1))*(ul1-geom[1])*(geom[1]-ll1)*(6.0/(ul2-ll2)/(ul2-ll2))*(ul2-geom[2])*(geom[2]-ll2);
                // faces perpendicular to Y-axis
                else if(side == 2 || side == 3)
                  specVal = DirichletData[aa][6]*(6.0/(ul1-ll1)/(ul1-ll1))*(ul1-geom[0])*(geom[0]-ll1)*(6.0/(ul2-ll2)/(ul2-ll2))*(ul2-geom[2])*(geom[2]-ll2);
                // faces perpendicular to Z-axis
                else
                  specVal = DirichletData[aa][6]*(6.0/(ul1-ll1)/(ul1-ll1))*(ul1-geom[0])*(geom[0]-ll1)*(6.0/(ul2-ll2)/(ul2-ll2))*(ul2-geom[1])*(geom[1]-ll2);
            }

            // multiply the velocity value with the time factor
            specVal *= timeFunction[0].prop;


              res = specVal - computeValue(dir, N);

              grad(0,0) = computeValue(0, dN_dx);
              grad(0,1) = computeValue(0, dN_dy);
              grad(0,2) = computeValue(0, dN_dz);

              grad(1,0) = computeValue(1, dN_dx);
              grad(1,1) = computeValue(1, dN_dy);
              grad(1,2) = computeValue(1, dN_dz);

              grad(2,0) = computeValue(2, dN_dx);
              grad(2,1) = computeValue(2, dN_dy);
              grad(2,2) = computeValue(2, dN_dz);

              pres      = computeValue(3, N);
              
              stress = mu*grad;
              stress(0,0) -= pres;
              stress(1,1) -= pres;
              stress(2,2) -= pres;

              trac = stress*normal ;

              for(ii=0;ii<totnlbf2;ii++)
              {
                fact = (dvol*PENALTY)* N[ii] ;

                TI = ndof*ii+dir;

                Flocal(TI) += fact*res;
                //fact *= af;

                for(jj=0;jj<totnlbf2;jj++)
                {
                  Klocal(TI, ndof*jj+dir) += fact * N[jj];
                }
              }

              // terms corresponding to Nitsche method for Stokes and Navier-Stokes
              if(isNitsche)
              {
                if(dir == 0)
                {
                  for(ii=0;ii<totnlbf2;ii++)
                  {
                    TI   = ndof*ii;
                    TIp1 = TI+1;
                    TIp2 = TI+2;
                    TIp3 = TI+3;

                    Ta = mu*(normal[0]*dN_dx[ii] + normal[1]*dN_dy[ii] + normal[2]*dN_dz[ii])*dvol;
                    bb1 = N[ii]*dvol;

                    for(jj=0;jj<totnlbf2;jj++)
                    {
                      TJ   = ndof*jj;
                      TJp1 = TJ+1;
                      TJp2 = TJ+2;
                      TJp3 = TJ+3;

                      Tb = af*mu*( normal[0]*dN_dx[jj] + normal[1]*dN_dy[jj] + normal[2]*dN_dz[jj] );

                      fact = bb1*(-normal[0]*N[jj]);

                      Klocal(TI, TJ)   -= (bb1*Tb);
                      Klocal(TI, TJp3) -= fact;

                      // Nitsche terms
                      Klocal(TI, TJ)   -= (Ta*af*N[jj])*NitscheFact;
                      Klocal(TIp3, TJ) -= fact*af*NitscheFact;
                    }

                    Flocal(TI)   -= (-bb1*trac[0]);
                    Flocal(TIp1) -= (0.0);
                    Flocal(TIp2) -= (0.0);

                    Flocal(TI)   -= (Ta*res)*NitscheFact;
                    Flocal(TIp3) -= (bb1*(-normal[0]*res))*NitscheFact;
                  }
                }
                if(dir == 1)
                {
                  for(ii=0;ii<totnlbf2;ii++)
                  {
                    TI   = ndof*ii;
                    TIp1 = TI+1;
                    TIp2 = TI+2;
                    TIp3 = TI+3;

                    Ta = mu*(normal[0]*dN_dx[ii] + normal[1]*dN_dy[ii] + normal[2]*dN_dz[ii])*dvol;
                    bb1 = N[ii]*dvol;

                    for(jj=0;jj<totnlbf2;jj++)
                    {
                      TJ   = ndof*jj;
                      TJp1 = TJ+1;
                      TJp2 = TJ+2;
                      TJp3 = TJ+3;

                      Tb = af*mu*( normal[0]*dN_dx[jj] + normal[1]*dN_dy[jj] + normal[2]*dN_dz[jj] );

                      fact = bb1*(-normal[1]*N[jj]);

                      Klocal(TIp1, TJp1)  -= (bb1*Tb);
                      Klocal(TIp1, TJp3)  -= fact;

                      // Nitsche terms
                      Klocal(TIp1, TJp1)  -= (Ta*af*N[jj])*NitscheFact;
                      Klocal(TIp3, TJp1)  -= fact*af*NitscheFact;
                    }

                    Flocal(TI)   -= (0.0);
                    Flocal(TIp1) -= (-bb1*trac[1]);
                    Flocal(TIp2) -= (0.0);

                    Flocal(TIp1) -= (Ta*res)*NitscheFact;
                    Flocal(TIp3) -= (bb1*(-normal[1]*res))*NitscheFact;
                  }
                }
                if(dir == 2)
                {
                  for(ii=0;ii<totnlbf2;ii++)
                  {
                    TI   = ndof*ii;
                    TIp1 = TI+1;
                    TIp2 = TI+2;
                    TIp3 = TI+3;

                    Ta = mu*(normal[0]*dN_dx[ii] + normal[1]*dN_dy[ii] + normal[2]*dN_dz[ii])*dvol;
                    bb1 = N[ii]*dvol;

                    for(jj=0;jj<totnlbf2;jj++)
                    {
                      TJ   = ndof*jj;
                      TJp1 = TJ+1;
                      TJp2 = TJ+2;
                      TJp3 = TJ+3;

                      Tb = af*mu*( normal[0]*dN_dx[jj] + normal[1]*dN_dy[jj] + normal[2]*dN_dz[jj] );

                      fact = bb1*(-normal[2]*N[jj]);

                      Klocal(TIp2, TJp2)  -= (bb1*Tb);
                      Klocal(TIp2, TJp3)  -= fact;

                      // Nitsche terms
                      Klocal(TIp2, TJp2)  -= (Ta*af*N[jj])*NitscheFact;
                      Klocal(TIp3, TJp2)  -= fact*af*NitscheFact;
                    }

                    Flocal(TI)   -= (0.0);
                    Flocal(TIp1) -= (0.0);
                    Flocal(TIp2) -= (-bb1*trac[2]);

                    Flocal(TIp2) -= (Ta*res)*NitscheFact;
                    Flocal(TIp3) -= (bb1*(-normal[2]*res))*NitscheFact;
                  }
                }
              } //if(isNitsche)
           //} // if(! GeomData->polyImm.within(xx, yy) )
        }// for(gp=0...
      } // for(aa=0;aa<DirichletData.size();aa++)
  } // if(DirichletData.size() > 0)

  return;
}




template<>
void TreeNode<3>::applyNeumannBCsCutFEMFluid(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
  if( NeumannData.size() > 0 )
  {
      int  aa, ii, jj, gp, nGauss, TI, TIp1, TIp2, index, dir, side;
      double  y0, y1, res, JacTemp;
      double  dvol, specVal, freq;

      VectorXd  N(totnlbf), NN(totnlbf);
      myPoint  param, normal, geom;
      double *gws;
      myPoint *gps;

      for(aa=0;aa<NeumannData.size();aa++)
      {
        side    = (int) (NeumannData[aa][0] - 1);
        dir     = (int) (NeumannData[aa][1] - 1);
        specVal = NeumannData[aa][2];

        normal = GeomData->boundaryNormals[side];

        if(domNums.size() > 1)
        {
          nGauss = BoundaryQuadrature[side].gausspoints.size() ;

          gps = &(BoundaryQuadrature[side].gausspoints[0]);
          gws = &(BoundaryQuadrature[side].gaussweights[0]);

          JacTemp = 1.0;
        }
        else
        {
          nGauss = GeomData->boundaryQuadrature3D[side].gausspoints.size();

          gps = &(GeomData->boundaryQuadrature3D[side].gausspoints[0]);
          gws = &(GeomData->boundaryQuadrature3D[side].gaussweights[0]);

          JacTemp = GeomData->boundaryJacobians[side][level];
        }

        for(gp=0; gp<nGauss; gp++)
        {
            param[0] = 0.5*(knotIncr[0] * gps[gp][0] + knotSum[0]);
            param[1] = 0.5*(knotIncr[1] * gps[gp][1] + knotSum[1]);
            param[2] = 0.5*(knotIncr[2] * gps[gp][2] + knotSum[2]);

            dvol = JacTemp * gws[gp] ;

            if(parent == NULL)
            {
              GeomData->computeBasisFunctions3D(knotBegin, knotIncr, param, N);
            }
            else
            {
              GeomData->computeBasisFunctions3D(knotBegin, knotIncr, param, NN);
              N = SubDivMat*NN;
            }

            geom[0] = GeomData->computeCoord(0, param[0]);
            geom[1] = GeomData->computeCoord(1, param[1]);
            geom[2] = GeomData->computeCoord(2, param[2]);

            //r = sqrt(xx*xx+yy*yy);
            //val = 1.0+log(2.0*r);
            //cout << xx << '\t' << yy << endl;

            specVal = NeumannData[aa][2];

            /*
            if(side == 110)
            {
              if(dir == 0)
              {
                if(abs(yy) >= 0.75)
                  specVal = 0.0;
                else
                  specVal = NeumannData[aa][2];
              }
            }
            */

            res = dvol * specVal * timeFunction[0].prop;
            //res = dvol * specVal * tanh(20.0*mpapTime.cur);
            //cout << specVal << '\t' << axsy << '\t' << dir << endl;

            //res *= (0.5*(1.0-cos(2.0*PI*SolnData->ElemProp.data[6]*mpapTime.cur)));
            //w1 = tanh(SolnData->ElemProp.data[5]*mpapTime.cur);

            freq = 2.0*PI*10.0;
            //res = dvol* specVal * (0.5*(1.0-cos(freq*mpapTime.cur)));
            //res = dvol * specVal * sin(freq*mpapTime.cur);
            //res *= tanh(50.0*timeFunction[0].prop);

            //cout << " yy " << yy << '\t' << res << endl;
            for(ii=0;ii<totnlbf2;ii++)
              Flocal(ndof*ii+dir) += (res*N[ii]);

        }// for(gp=0...
      } // for(aa=0;aa<NeumannData.size();aa++)
  } // if(NeumannData.size() > 0)

  return;
}



template<>
void TreeNode<3>::applyDerivativeBCsCutFEMFluid(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
  if( NeumannData.size() > 0 )
  {
      int  aa, ii, jj, gp, nGauss, TI, TIp1, TIp2, index, dir, side;
      double  y0, y1, res, JacTemp;
      double  dvol, specVal, freq;

      VectorXd  N(totnlbf), NN(totnlbf);
      myPoint  param, normal, geom;
      double *gws;
      myPoint *gps;

      for(aa=0;aa<NeumannData.size();aa++)
      {
        side    = (int) (NeumannData[aa][0] - 1);
        dir     = (int) (NeumannData[aa][1] - 1);
        specVal = NeumannData[aa][2];

        normal = GeomData->boundaryNormals[side];

        if(domNums.size() > 1)
        {
          nGauss = BoundaryQuadrature[side].gausspoints.size() ;

          gps = &(BoundaryQuadrature[side].gausspoints[0]);
          gws = &(BoundaryQuadrature[side].gaussweights[0]);

          JacTemp = 1.0;
        }
        else
        {
          nGauss = GeomData->boundaryQuadrature3D[side].gausspoints.size();

          gps = &(GeomData->boundaryQuadrature3D[side].gausspoints[0]);
          gws = &(GeomData->boundaryQuadrature3D[side].gaussweights[0]);

          JacTemp = GeomData->boundaryJacobians[side][level];
        }

        for(gp=0; gp<nGauss; gp++)
        {
            param[0] = 0.5*(knotIncr[0] * gps[gp][0] + knotSum[0]);
            param[1] = 0.5*(knotIncr[1] * gps[gp][1] + knotSum[1]);
            param[2] = 0.5*(knotIncr[2] * gps[gp][2] + knotSum[2]);

            dvol = JacTemp * gws[gp] ;

            if(parent == NULL)
            {
              GeomData->computeBasisFunctions3D(knotBegin, knotIncr, param, N);
            }
            else
            {
              GeomData->computeBasisFunctions3D(knotBegin, knotIncr, param, NN);
              N = SubDivMat*NN;
            }

            geom[0] = GeomData->computeCoord(0, param[0]);
            geom[1] = GeomData->computeCoord(1, param[1]);
            geom[2] = GeomData->computeCoord(2, param[2]);

            //r = sqrt(xx*xx+yy*yy);
            //val = 1.0+log(2.0*r);
            //cout << xx << '\t' << yy << endl;

            specVal = NeumannData[aa][2];

            /*
            if(side == 110)
            {
              if(dir == 0)
              {
                if(abs(yy) >= 0.75)
                  specVal = 0.0;
                else
                  specVal = NeumannData[aa][2];
              }
            }
            */

            res = dvol * specVal * timeFunction[0].prop;
            //res = dvol * specVal * tanh(20.0*mpapTime.cur);
            //cout << specVal << '\t' << axsy << '\t' << dir << endl;

            //res *= (0.5*(1.0-cos(2.0*PI*SolnData->ElemProp.data[6]*mpapTime.cur)));
            //w1 = tanh(SolnData->ElemProp.data[5]*mpapTime.cur);

            freq = 2.0*PI*10.0;
            //res = dvol* specVal * (0.5*(1.0-cos(freq*mpapTime.cur)));
            //res = dvol * specVal * sin(freq*mpapTime.cur);
            //res *= tanh(50.0*timeFunction[0].prop);

            //cout << " yy " << yy << '\t' << res << endl;
            for(ii=0;ii<totnlbf2;ii++)
              Flocal(ndof*ii+dir) += (res*N[ii]);

        }// for(gp=0...
      } // for(aa=0;aa<NeumannData.size();aa++)
  } // if(NeumannData.size() > 0)

  return;
}








