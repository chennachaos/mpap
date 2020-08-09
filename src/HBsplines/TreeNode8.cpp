
#include "TreeNode.h"
#include "MpapTime.h"
#include "TimeFunction.h"
#include "Functions.h"
#include "SolutionData.h"
#include "BasisFunctionsBSpline.h"
#include "myDataIntegrateCutFEM.h"
#include "myPoly.h"
#include "stabilisationRoutines.h"
#include "myLine.h"
#include "myTria.h"


extern  MpapTime  mpapTime;
extern List<TimeFunction> timeFunction;

using namespace myGeom;



/*
    fact = sqrt(volume/(hx*hy));

    hx *= fact;
    hy *= fact;

    matJ.setZero();
    matJ(0,0) = hx*0.5;
    matJ(1,1) = hy*0.5;

    matJinv = matJ.inverse();

    matG.setZero();
    matG(0,0) = matJinv(0,0);
    matG(0,1) = matJinv(0,1);
    matG(1,0) = matJinv(1,0);
    matG(1,1) = matJinv(1,1);

    matG = matG.transpose() * matG;
*/

template<>
void TreeNode<1>::calcStiffnessAndResidualCutFEMFluid(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
  return;
}



/*
template<>
void TreeNode<2>::calcStiffnessAndResidualCutFEMFluid(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
    // GFEM for Stokes
    // with proper stabilisation
    // with Newton-Raphson
    ///////////////////////////////////////

    int ii, jj, gp, nGauss, tempId, count=0;
    int TI, TIp1, TIp2, TJ, TJp1, TJp2;

    double  JacTemp, Jac, dvol, stabParam, CI=4.0;
    double  fact, fact2, b1, b2, b3, b4, b5, b6, b7, b8, acceFact, dt, bforce[2];
    double  pres, Da, Db, af, am, muTaf, rad, urdr, urdr2, h2, h, tau[3];

    double  beta[6]; get_stabilisation_beta_wulf(beta);

    VectorXd  NN(totnlbf), dNN_dx(totnlbf), d2NN_dx2(totnlbf), dNN_dy(totnlbf), d2NN_dy2(totnlbf);
    VectorXd  N, dN_dx, d2N_dx2, dN_dy, d2N_dy2, d2N, velTemp(3);
    VectorXd  res(3), res2(2), dp(2), Du(2), vel(2), velDot(2), force(2), gradTvel(2), rStab(3);
    MatrixXd  Dj(2, 3), grad(2,2), gradN(2,2), stress(2,2);
    MatrixXd  matB(2,2), matBinv(2,2), matG(3,3);
    myPoint  param, geom, velPrev;
    Dj.setZero();


    bool   axsy = ((int)elmDat[2] == 1);
    double  rho = elmDat[3];
    double  mu  = elmDat[4];
    bforce[0]   = elmDat[5];
    bforce[1]   = elmDat[6];


    af = SolnData->td(2);
    am = SolnData->td(1);
    acceFact = am*SolnData->td(9);
    dt = mpapTime.dt;

    muTaf = mu*af;

    double *gws;
    myPoint *gps;

    double  hx = bbox.maxBB[0]-bbox.minBB[0];
    double  hy = bbox.maxBB[1]-bbox.minBB[1];

    volume = hx*hy;

    if(domNums.size() > 1)
    {
      nGauss = Quadrature.gausspoints.size();
      
      gps = &(Quadrature.gausspoints[0]);
      gws = &(Quadrature.gaussweights[0]);
      
      JacTemp = 1.0;

      volume = 0.0;
      for(gp=0; gp<nGauss; gp++)
        volume += gws[gp];
    }
    else
    {
      nGauss = GeomData->gausspoints.size();

      gps = &(GeomData->gausspoints[0]);
      gws = &(GeomData->gaussweights[0]);

      JacTemp = JacMultElem;
    }

    //cout << " nGauss " << nGauss << '\t' << tempId << endl;

    fact = volume/(hx*hy);

    fact = sqrt(fact);

    hx *= fact;
    hy *= fact;
    
    matB.setZero();
    matB(0,0) = hx*0.5;
    matB(1,1) = hy*0.5;

    matBinv = matB.inverse();

    matG.setZero();
    matG(0,0) = matBinv(0,0);
    matG(0,1) = matBinv(0,1);
    matG(1,0) = matBinv(1,0);
    matG(1,1) = matBinv(1,1);

    matG = matG.transpose() * matG;

//    cout << volume << '\t' << volume << endl;
    h2 = 4.0*volume/PI;
    h = sqrt(h2);

    stabParam = h2/(4.0*mu);
    //stabParam /= degree[0]/degree[1];
    ////stabParam /= rho;
    //stabParam *= rho;

    tau[0] = elmDat[8]*stabParam;  // SUPG
    tau[1] = elmDat[9]*stabParam;  // PSPG
    tau[2] = elmDat[10]*stabParam; // LSIC

    //cout << tau[0] << '\t' << tau[1] << '\t' << tau[2] << endl;

    //KimMoinFlow  analy(rho, mu);
    //Kovasznay  analy;
    //Stokes2DEx1  analy;

    count=0;
    for(gp=0; gp<nGauss; gp++)
    {
        param[0]  = 0.5*(knotIncr[0] * gps[gp][0] + knotSum[0]);
        param[1]  = 0.5*(knotIncr[1] * gps[gp][1] + knotSum[1]);

        dvol = gws[gp] * JacTemp;

        geom[0] = GeomData->computeCoord(0, param[0]);
        geom[1] = GeomData->computeCoord(1, param[1]);

          GeomData->computeBasisFunctions2D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy, d2NN_dx2, d2NN_dy2);

          if(parent == NULL)
          {
            //cout << " parent is NULL " << endl;
            N = NN;
            dN_dx = dNN_dx;
            dN_dy = dNN_dy;
            d2N_dx2 = d2NN_dx2;
            d2N_dy2 = d2NN_dy2;
          }
          else
          {
            //cout << " parent is not NULL " << endl;
            N = SubDivMat*NN;
            dN_dx = SubDivMat*dNN_dx;
            dN_dy = SubDivMat*dNN_dy;
            d2N_dx2 = SubDivMat*d2NN_dx2;
            d2N_dy2 = SubDivMat*d2NN_dy2;
          }

        //cout <<  " gpDomainId = " << gpDomainId << '\t' << totnlbf2 << endl;
        //printVector(GlobalBasisFuncs);

        d2N = d2N_dx2 + d2N_dy2;

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

            b1     = SolnData->var1Prev(TI);
            b2     = SolnData->var1Prev(TIp1);

            velPrev(0) += ( b1 * N(ii) );
            velPrev(1) += ( b2 * N(ii) );

            b1     = SolnData->var1Cur(TI);
            b2     = SolnData->var1Cur(TIp1);
            b3     = SolnData->var1(TIp2);

            vel(0) += ( b1 * N(ii) );
            vel(1) += ( b2 * N(ii) );

            grad(0,0) += ( b1 * dN_dx(ii) );
            grad(0,1) += ( b1 * dN_dy(ii) );
            grad(1,0) += ( b2 * dN_dx(ii) );
            grad(1,1) += ( b2 * dN_dy(ii) );

            Du(0)  += ( b1 * d2N(ii) );
            Du(1)  += ( b2 * d2N(ii) );

            pres   += ( b3 * N(ii) );
            dp(0)  += ( b3 * dN_dx(ii) );
            dp(1)  += ( b3 * dN_dy(ii) );

            b1     = SolnData->var1DotCur(TI);
            b2     = SolnData->var1DotCur(TIp1);

            velDot(0) += ( b1 * N(ii) );
            velDot(1) += ( b2 * N(ii) );
          }

          // this is pseudo-stress
          //stress = mu*(grad+grad.transpose());
          stress = mu*grad;
          stress(0,0) -= pres;
          stress(1,1) -= pres;

          force.setZero();

          //force(0) = analy.computeForce(0, geom[0], geom[1]);
          //force(1) = analy.computeForce(1, geom[0], geom[1]);

          //force(0) = GeomData->analyDBC->computeForce(0, xx, yy);
          //force(1) = GeomData->analyDBC->computeForce(0, xx, yy);
          //force = 1.0;
          //force = 0.0;
          //cout << force(0) << '\t' << force(1) << endl;

          //force[0] -= bforce[0];
          //force[1] -= bforce[1];

          res2(0) = rho*(velDot(0) - force(0)) ;
          res2(1) = rho*(velDot(1) - force(1)) ;

          rStab(0) = res2(0) - mu*Du(0) + dp(0) ;
          rStab(1) = res2(1) - mu*Du(1) + dp(1) ;

          if(axsy)
          {
            rad = geom[0];

            urdr  = vel(0)/rad;
            urdr2 = urdr/rad;

            dvol *= (2.0*PI*rad);
            
            rStab(0) -= mu*(grad(0,0)/rad - urdr2 );
            rStab(1) -= mu*(grad(1,0)/rad );
          }

          // evaluate stabilisation parameters
          //
          velTemp(0) = velPrev(0);
          velTemp(1) = velPrev(1);
          velTemp(2) = 0.0;

          //evaluateStabParams_algo1(&velTemp(0), h, rho, mu, dt,  beta, tau);

          //evaluateStabParams_algo2(&velTemp(0), h, rho, mu, dt,  beta, tau);

          //evaluateStabParams_algo3(velTemp, matG, dt, rho, mu, CI, tau);

          //if( abs(mpapTime.cur - mpapTime.dt) < 1.0e-10 )
            //tau[0] = 0.0;
          //else
            //tau[0] *= elmDat[8];  // SUPG

          //tau[1] *= elmDat[9];  // PSPG
          //tau[2] *= elmDat[10]; // LSIC

          for(ii=0;ii<totnlbf2;ii++)
          {
            TI   = ndof*ii;
            TIp1 = TI+1;
            TIp2 = TI+2;

            b1 = dN_dx[ii]*dvol;
            b2 = dN_dy[ii]*dvol;
            b4 = N[ii]*dvol;

            b5 = muTaf*b1;
            b6 = muTaf*b2;
            b8 = af*b4;

            for(jj=0;jj<totnlbf2;jj++)
            {
              TJ   = ndof*jj;
              TJp1 = TJ+1;
              TJp2 = TJ+2;

              fact2 = rho*acceFact*N(jj);

              // time acceleration term
              fact = b4*fact2 ;

              // diffusion term
              fact += ( b5*dN_dx(jj)+b6*dN_dy(jj) );

              Klocal(TI,   TJ)   += fact;
              Klocal(TIp1, TJp1) += fact;

              //Klocal(TI,   TJ)   += b5*dN_dx(jj);
              //Klocal(TI,   TJp1) += b6*dN_dx(jj);
              //Klocal(TIp1, TJ)   += b5*dN_dy(jj);
              //Klocal(TIp1, TJp1) += b6*dN_dy(jj);

              // pressure term
              Klocal(TI,   TJp2) -= (b1*N(jj));
              Klocal(TIp1, TJp2) -= (b2*N(jj));

              // continuity equation
              Klocal(TIp2, TJ)   += (b8*dN_dx(jj));
              Klocal(TIp2, TJp1) += (b8*dN_dy(jj));

              // SUPG and PSPG stabilisation terms
              fact2 -= ( muTaf*d2N(jj) );

              Dj(0,0) = fact2;
              Dj(0,1) = 0.0;
              Dj(0,2) = dN_dx(jj);
              Dj(1,0) = 0.0;
              Dj(1,1) = fact2;
              Dj(1,2) = dN_dy(jj);
              
              if(axsy)
              {
                Dj(0,0) -= muTaf*(dN_dx(jj)/rad - N(jj)/rad/rad);
                Dj(1,1) -= muTaf*(dN_dx(jj)/rad);
              }

              // PSPG
              Klocal(TIp2, TJ)   += (b1*Dj(0,0) + b2*Dj(1,0))*tau[1];
              Klocal(TIp2, TJp1) += (b1*Dj(0,1) + b2*Dj(1,1))*tau[1];
              Klocal(TIp2, TJp2) += (b1*Dj(0,2) + b2*Dj(1,2))*tau[1];

              // LSIC stabilisation

              fact2 = rho*af*tau[2];
              Klocal(TI,   TJ)   += (b1*dN_dx(jj))*fact2;
              Klocal(TI,   TJp1) += (b1*dN_dy(jj))*fact2;

              Klocal(TIp1, TJ)   += (b2*dN_dx(jj))*fact2;
              Klocal(TIp1, TJp1) += (b2*dN_dy(jj))*fact2;

              if(axsy)
              {
                  // diffusion term
                  Klocal(TI, TJ)     += (b4 * (mu/rad/rad) * (af*N(jj)) );
                  Klocal(TI, TJp2)   -= (b4 * N(jj)/rad);

                  // continuity equation
                  Klocal(TIp2, TJ)   += (b4 * af*N(jj)/rad);
              }
            }

            Flocal(TI)   -= (b4*res2(0) + b1*stress(0,0) + b2*stress(0,1) );
            Flocal(TIp1) -= (b4*res2(1) + b1*stress(1,0) + b2*stress(1,1) );
            Flocal(TIp2) -= (b4*grad.trace());

            // PSPG stabilisation terms
            Flocal(TIp2) -= (tau[1]*(b1*rStab(0)+b2*rStab(1)));

            // LSIC stabilisation terms
            
            fact2 = tau[2]*rho*grad.trace();
            
            Flocal(TI)   -= b1*fact2;
            Flocal(TIp1) -= b2*fact2;

            if(axsy)
            {
                Flocal(TI)   -= (b4 * (mu/rad/rad) * vel(0) );
                Flocal(TI)   += (b4 * pres/rad);
                Flocal(TIp2) -= (b4 * vel(0)/rad);
            }
          } // for(ii=0;ii<totnlbf2;ii++)
      //} //if(within)
    }//gp

    return;
}
*/




/*
template<>
void TreeNode<2>::calcStiffnessAndResidualCutFEMFluid(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
    // GFEM for Stokes
    // with proper stabilisation
    // without Newton-Raphson
    ///////////////////////////////////////

    int ii, jj, gp, nGauss, tempId, count=0;
    int TI, TIp1, TIp2, TJ, TJp1, TJp2;

    double  JacTemp, Jac, dvol, stabParam, CI=4.0;
    double  fact, fact2, b1, b2, b3, b4, b5, b6, b7, b8, acceFact, dt, bforce[2];
    double  pres, Da, Db, af, am, muTaf, rad, urdr, urdr2, h2, h, tau[3];

    double  beta[6]; get_stabilisation_beta_wulf(beta);

    VectorXd  NN(totnlbf), dNN_dx(totnlbf), d2NN_dx2(totnlbf), dNN_dy(totnlbf), d2NN_dy2(totnlbf);
    VectorXd  N, dN_dx, d2N_dx2, dN_dy, d2N_dy2, d2N, velTemp(3);
    VectorXd  res(3), res2(2), dp(2), Du(2), vel(2), velDot(2), force(2), gradTvel(2), rStab(3);
    MatrixXd  Dj(2, 3), grad(2,2), gradN(2,2), stress(2,2), matG(3,3);
    myPoint  param, geom, velPrev;
    Dj.setZero();


    bool   axsy = ((int)elmDat[2] == 1);
    double  rho = elmDat[3];
    double  mu  = elmDat[4];
    bforce[0]   = elmDat[5];
    bforce[1]   = elmDat[6];


    af = SolnData->td(2);
    am = SolnData->td(1);
    acceFact = am*SolnData->td(9);
    dt = mpapTime.dt;

    muTaf = mu*af;

    double *gws;
    myPoint *gps;

    double  hx = bbox.maxBB[0]-bbox.minBB[0];
    double  hy = bbox.maxBB[1]-bbox.minBB[1];

    volume = hx*hy;

    matG.setZero();
    matG(0,0) = 4.0/hx/hx;
    matG(1,1) = 4.0/hy/hy;

    if(domNums.size() > 1)
    {
      nGauss = Quadrature.gausspoints.size();
      
      gps = &(Quadrature.gausspoints[0]);
      gws = &(Quadrature.gaussweights[0]);
      
      JacTemp = 1.0;

      volume = 0.0;
      for(gp=0; gp<nGauss; gp++)
        volume += gws[gp];

      // For 2D problem
      // fact = sqrt(Vc/V); and fact = fact*fact;  ---->  fact = Vc/V;
      // For 3D problem
      // fact = (Vc/V)^(1/3). So, do fact = fact*fact;

      fact = volume/(hx*hy);

      matG(0,0) /= fact;
      matG(1,1) /= fact;
    }
    else
    {
      nGauss = GeomData->gausspoints.size();

      gps = &(GeomData->gausspoints[0]);
      gws = &(GeomData->gaussweights[0]);
      
      JacTemp = JacMultElem;
    }

    //cout << volume << '\t' << volume << endl;
    //h2 = 4.0*volume/PI;
    //h = sqrt(h2);
    //stabParam = h2/(4.0*mu);
    //stabParam /= degree[0]/degree[1];
    ////stabParam /= rho;
    //stabParam *= rho;

    //tau[0] = elmDat[8]*stabParam;  // SUPG
    //tau[1] = elmDat[9]*stabParam;  // PSPG
    //tau[2] = elmDat[10]*stabParam; // LSIC

    //KimMoinFlow  analy(rho, mu);
    //Kovasznay  analy;

    count=0;
    for(gp=0; gp<nGauss; gp++)
    {
        param[0]  = 0.5*(knotIncr[0] * gps[gp][0] + knotSum[0]);
        param[1]  = 0.5*(knotIncr[1] * gps[gp][1] + knotSum[1]);

        dvol = gws[gp] * JacTemp;

        geom[0] = GeomData->computeCoord(0, param[0]);
        geom[1] = GeomData->computeCoord(1, param[1]);

        //cout << uu << '\t' << vv << endl;
        //cout << xx << '\t' << yy << endl;

        //cout << gp << '\t' << gps[gp][0] << '\t' << gps[gp][1] << '\t' << gws[gp] << '\t' << dvol << endl;

          GeomData->computeBasisFunctions2D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy, d2NN_dx2, d2NN_dy2);

          if(parent == NULL)
          {
            //cout << " parent is NULL " << endl;
            N = NN;
            dN_dx = dNN_dx;
            dN_dy = dNN_dy;
            d2N_dx2 = d2NN_dx2;
            d2N_dy2 = d2NN_dy2;
          }
          else
          {
            //cout << " parent is not NULL " << endl;
            N = SubDivMat*NN;
            dN_dx = SubDivMat*dNN_dx;
            dN_dy = SubDivMat*dNN_dy;
            d2N_dx2 = SubDivMat*d2NN_dx2;
            d2N_dy2 = SubDivMat*d2NN_dy2;
          }

        //cout <<  " gpDomainId = " << gpDomainId << '\t' << totnlbf2 << endl;
        //printVector(GlobalBasisFuncs);

        d2N = d2N_dx2 + d2N_dy2;

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

            b1     = SolnData->var1Prev(TI);
            b2     = SolnData->var1Prev(TIp1);

            velPrev(0) += ( b1 * N(ii) );
            velPrev(1) += ( b2 * N(ii) );

            b1     = SolnData->var1Cur(TI);
            b2     = SolnData->var1Cur(TIp1);
            b3     = SolnData->var1(TIp2);

            vel(0) += ( b1 * N(ii) );
            vel(1) += ( b2 * N(ii) );

            grad(0,0) += ( b1 * dN_dx(ii) );
            grad(0,1) += ( b1 * dN_dy(ii) );
            grad(1,0) += ( b2 * dN_dx(ii) );
            grad(1,1) += ( b2 * dN_dy(ii) );

            Du(0)  += ( b1 * d2N(ii) );
            Du(1)  += ( b2 * d2N(ii) );

            pres   += ( b3 * N(ii) );
            dp(0)  += ( b3 * dN_dx(ii) );
            dp(1)  += ( b3 * dN_dy(ii) );

            b1     = SolnData->var1DotCur(TI);
            b2     = SolnData->var1DotCur(TIp1);

            velDot(0) += ( b1 * N(ii) );
            velDot(1) += ( b2 * N(ii) );
          }

          // this is pseudo-stress
          //stress = mu*(grad+grad.transpose());
          stress = mu*grad;
          stress(0,0) -= pres;
          stress(1,1) -= pres;

          force.setZero();

          res2(0) = rho*(velPrev(0)/dt ) ;
          res2(1) = rho*(velPrev(1)/dt ) ;

          rStab(0) = res2(0) ;
          rStab(1) = res2(1) ;

          // evaluate stabilisation parameters
          //
          velTemp(0) = velPrev(0);
          velTemp(1) = velPrev(1);
          velTemp(2) = 0.0;

          //evaluateStabParams_algo1(&velTemp(0), h, rho, mu, dt,  beta, tau);

          //evaluateStabParams_algo2(&velTemp(0), h, rho, mu, dt,  beta, tau);

          evaluateStabParams_algo3(velTemp, matG, dt, rho, mu, CI, tau);

          //if( abs(mpapTime.cur - mpapTime.dt) < 1.0e-10 )
            //tau[0] = 0.0;
          //else
            tau[0] *= elmDat[8];  // SUPG

          tau[1] *= elmDat[9];  // PSPG
          tau[2] *= elmDat[10]; // LSIC

          for(ii=0;ii<totnlbf2;ii++)
          {
            TI   = ndof*ii;
            TIp1 = TI+1;
            TIp2 = TI+2;

            b1 = dN_dx[ii]*dvol;
            b2 = dN_dy[ii]*dvol;
            b4 = N[ii]*dvol;

            b5 = muTaf*b1;
            b6 = muTaf*b2;
            b8 = af*b4;

            for(jj=0;jj<totnlbf2;jj++)
            {
              TJ   = ndof*jj;
              TJp1 = TJ+1;
              TJp2 = TJ+2;

              fact2 = rho*acceFact*N(jj);

              // time acceleration term
              fact = b4*fact2 ;

              // diffusion term
              fact += ( b5*dN_dx(jj)+b6*dN_dy(jj) );

              Klocal(TI,   TJ)   += fact;
              Klocal(TIp1, TJp1) += fact;

              // pressure term
              Klocal(TI,   TJp2) -= (b1*N(jj));
              Klocal(TIp1, TJp2) -= (b2*N(jj));

              // continuity equation
              Klocal(TIp2, TJ)   += (b8*dN_dx(jj));
              Klocal(TIp2, TJp1) += (b8*dN_dy(jj));

              // SUPG and PSPG stabilisation terms
              fact2 -= ( muTaf*d2N(jj) );

              Dj(0,0) = fact2;
              Dj(0,1) = 0.0;
              Dj(0,2) = dN_dx(jj);
              Dj(1,0) = 0.0;
              Dj(1,1) = fact2;
              Dj(1,2) = dN_dy(jj);

              // PSPG
              Klocal(TIp2, TJ)   += (b1*Dj(0,0) + b2*Dj(1,0))*tau[1];
              Klocal(TIp2, TJp1) += (b1*Dj(0,1) + b2*Dj(1,1))*tau[1];
              Klocal(TIp2, TJp2) += (b1*Dj(0,2) + b2*Dj(1,2))*tau[1];

            }

            Flocal(TI)   += (b4*res2(0)  );
            Flocal(TIp1) += (b4*res2(1)  );

            // PSPG stabilisation terms
            Flocal(TIp2) += (tau[1]*(b1*rStab(0)+b2*rStab(1)));

          } // for(ii=0;ii<totnlbf2;ii++)
      //} //if(within)
    }//gp

    return;
}
*/



/*
template<>
void TreeNode<2>::calcStiffnessAndResidualCutFEMFluid(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
    // basis functions are computed based not computed everytime

    // GFEM for Navier-Stokes
    // with proper stabilisation
    // fully-implicit formulation
    // 
    ///////////////////////////////////////

    int  ii=0, jj=0, gp=0, nGauss=0, tempId=0, count=0;
    int  TI=0, TIp1=0, TIp2=0, TJ=0, TJp1=0, TJp2=0;

    double  JacTemp=0.0, Jac=0.0, dvol=0.0, stabParam=0.0, CI=7.0;
    double  fact=0.0, fact2=0.0, b1=0.0, b2=0.0, b3=0.0, b4=0.0, b5=0.0, b6=0.0, b7=0.0, b8=0.0;
    double  pres=0.0, Da=0.0, Db=0.0, rad=0.0, urdr=0.0, urdr2=0.0, h2=0.0, h=0.0, tau[3];

    //double beta[6]; get_stabilisation_beta_wulf(beta);


    VectorXd  NN(totnlbf), dNN_dx(totnlbf), d2NN_dx2(totnlbf), dNN_dy(totnlbf), d2NN_dy2(totnlbf);
    VectorXd  N(totnlbf), dN_dx(totnlbf), d2N_dx2(totnlbf), dN_dy(totnlbf), d2N_dy2(totnlbf), d2N(totnlbf);
    VectorXd  res(3), res2(2), dp(2), Du(2), vel(2), velDot(2), velTemp(3);
    VectorXd  force(2), gradTvel(2), rStab(3), velPrev(2);
    MatrixXd  Dj(2, 3), grad(2,2), gradN(2,2), stress(2,2), gradPrev(2,2);
    MatrixXd  matJ(2,2), matJinv(2,2), matG(3,3);
    myPoint  param, geom;
    Dj.setZero();
    force.setZero();


    bool   axsy = ((int)elmDat[2] == 1);
    double  rho = elmDat[3];
    double  mu  = elmDat[4];
    double  bforce[2] = {elmDat[5], elmDat[6]};
    double  af = SolnData->td(2);
    double  am = SolnData->td(1);
    double  acceFact = am*SolnData->td(9);
    double  dt = mpapTime.dt;
    double  muTaf = mu*af;

    double *gws;
    myPoint *gps;

    double  hx = bbox.maxBB[0]-bbox.minBB[0];
    double  hy = bbox.maxBB[1]-bbox.minBB[1];

    double volume = hx*hy;

    matG.setZero();
    matG(0,0) = 4.0/hx/hx;
    matG(1,1) = 4.0/hy/hy;


    if(domNums.size() > 1) // cut cell
    {
      nGauss = Quadrature.gausspoints.size();

      gps = &(Quadrature.gausspoints[0]);
      gws = &(Quadrature.gaussweights[0]);

      JacTemp = 1.0;

      volume = 0.0;
      for(gp=0; gp<nGauss; gp++)
        volume += gws[gp];

      // For 2D problem
      // fact = sqrt(Vc/V); and fact = fact*fact;  ---->  fact = Vc/V;
      // For 3D problem
      // fact = (Vc/V)^(1/3). So, do fact = fact*fact;

      fact = volume/(hx*hy);

      matG(0,0) /= fact;
      matG(1,1) /= fact;
    }
    else
    {
      nGauss = GeomData->gausspoints.size();

      gps = &(GeomData->gausspoints[0]);
      gws = &(GeomData->gaussweights[0]);

      JacTemp = JacMultElem;
    }

    //double h2 = 4.0*volume/PI;
    //double h = sqrt(h2);
    //double stabParam = h2/(4.0*mu);
    //stabParam /= degree[0]/degree[1];
    ////stabParam /= rho;
    //stabParam *= rho;

    //tau[0] = elmDat[8]*stabParam;  // SUPG
    //tau[1] = elmDat[9]*stabParam;  // PSPG
    //tau[2] = elmDat[10]*stabParam; // LSIC

    //cout << tau[0] << '\t' << tau[1] << '\t' << tau[2] << endl;

    //KimMoinFlow  analy(rho, mu);
    //Kovasznay  analy;

    count=0;
    for(gp=0; gp<nGauss; gp++)
    {
        param[0]  = 0.5*(knotIncr[0] * gps[gp][0] + knotSum[0]);
        param[1]  = 0.5*(knotIncr[1] * gps[gp][1] + knotSum[1]);

        dvol = gws[gp] * JacTemp;

        geom[0] = GeomData->computeCoord(0, param[0]); // radius for axsy problems
        geom[1] = GeomData->computeCoord(1, param[1]);

        if(parent == NULL)  //cout << " parent is NULL " << endl;
        {
          GeomData->computeBasisFunctions2D(knotBegin, knotIncr, param, N, dN_dx, dN_dy, d2N_dx2, d2N_dy2);
        }
        else
        {
          GeomData->computeBasisFunctions2D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy, d2NN_dx2, d2NN_dy2);
          //cout << " parent is not NULL " << endl;
          N = SubDivMat*NN;
          dN_dx = SubDivMat*dNN_dx;
          dN_dy = SubDivMat*dNN_dy;
          d2N_dx2 = SubDivMat*d2NN_dx2;
          d2N_dy2 = SubDivMat*d2NN_dy2;
        }

        d2N = d2N_dx2 + d2N_dy2;


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

            b1     = SolnData->var1Prev(TI);
            b2     = SolnData->var1Prev(TIp1);

            velPrev(0) += ( b1 * N(ii) );
            velPrev(1) += ( b2 * N(ii) );

            b1     = SolnData->var1Cur(TI);
            b2     = SolnData->var1Cur(TIp1);
            b3     = SolnData->var1(TIp2);

            vel(0) += ( b1 * N(ii) );
            vel(1) += ( b2 * N(ii) );

            grad(0,0) += ( b1 * dN_dx(ii) );
            grad(0,1) += ( b1 * dN_dy(ii) );
            grad(1,0) += ( b2 * dN_dx(ii) );
            grad(1,1) += ( b2 * dN_dy(ii) );

            Du(0)  += ( b1 * d2N(ii) );
            Du(1)  += ( b2 * d2N(ii) );

            pres   += ( b3 * N(ii) );
            dp(0)  += ( b3 * dN_dx(ii) );
            dp(1)  += ( b3 * dN_dy(ii) );

            b1     = SolnData->var1DotCur(TI);
            b2     = SolnData->var1DotCur(TIp1);

            velDot(0) += ( b1 * N(ii) );
            velDot(1) += ( b2 * N(ii) );
          }

          // this is pseudo-stress
          //stress = mu*(grad+grad.transpose());
          stress = mu*grad;
          stress(0,0) -= pres;
          stress(1,1) -= pres;

          force.setZero();

          //force(0) = analy.computeForce(0, geom[0], geom[1]);
          //force(1) = analy.computeForce(1, geom[0], geom[1]);

          //force(0) = GeomData->analyDBC->computeForce(0, xx, yy);
          //force(1) = GeomData->analyDBC->computeForce(0, xx, yy);
          //force = 1.0;
          //force = 0.0;
          //cout << force(0) << '\t' << force(1) << endl;

          //force[0] -= bforce[0];
          //force[1] -= bforce[1];

          gradTvel = grad*vel;

          res2(0) = rho*(velDot(0) + gradTvel(0) - force(0)) ;
          res2(1) = rho*(velDot(1) + gradTvel(1) - force(1)) ;

          rStab(0) = res2(0) - mu*Du(0) + dp(0) ;
          rStab(1) = res2(1) - mu*Du(1) + dp(1) ;

          if(axsy)
          {
            rad = geom[0];

            urdr  = vel(0)/rad;
            urdr2 = urdr/rad;

            dvol *= (2.0*PI*rad);

            rStab(0) -= mu*(grad(0,0)/rad - urdr2 );
            rStab(1) -= mu*(grad(1,0)/rad );
          }

          // evaluate stabilisation parameters
          //
          velTemp(0) = velPrev(0);
          velTemp(1) = velPrev(1);
          velTemp(2) = 0.0;

          //evaluateStabParams_algo1(&velTemp(0), h, rho, mu, dt,  beta, tau);

          //evaluateStabParams_algo2(&velTemp(0), h, rho, mu, dt,  beta, tau);

          evaluateStabParams_algo3(velTemp, matG, dt, rho, mu, CI, tau);

<<<<<<< HEAD:HBsplines/TreeNode9.cpp
          if( abs(mpapTime.cur - mpapTime.dt) < 1.0e-10 )
            tau[0] = 0.0;
          else
            tau[0] *= elmDat[8];  // SUPG
=======
          //if( abs(mpapTime.cur - mpapTime.dt) < 1.0e-10 )
            //tau[0] = 0.0;
          //else
            //tau[0] *= elmDat[8];  // SUPG
>>>>>>> collabchandan:src/HBsplines/TreeNode8.cpp

          tau[0] *= elmDat[8];  // SUPG
          tau[1] *= elmDat[9];  // PSPG
          tau[2] *= elmDat[10]; // LSIC

          for(ii=0;ii<totnlbf2;ii++)
          {
            TI   = ndof*ii;
            TIp1 = TI+1;
            TIp2 = TI+2;

            b1 = dN_dx[ii]*dvol;
            b2 = dN_dy[ii]*dvol;
            b4 = N[ii]*dvol;

            b5 = muTaf*b1;
            b6 = muTaf*b2;
            b8 = af*b4;

            Da = rho*(vel(0)*b1 + vel(1)*b2)*tau[0];

            for(jj=0;jj<totnlbf2;jj++)
            {
              TJ   = ndof*jj;
              TJp1 = TJ+1;
              TJp2 = TJ+2;

              fact2 = rho*acceFact*N(jj);

              // time acceleration term
              fact = b4*fact2 ;

              // diffusion term
              fact += ( b5*dN_dx(jj)+b6*dN_dy(jj) );

              Klocal(TI,   TJ)   += fact;
              Klocal(TIp1, TJp1) += fact;

              //Klocal(TI,   TJ)   += b5*dN_dx(jj);
              //Klocal(TI,   TJp1) += b6*dN_dx(jj);
              //Klocal(TIp1, TJ)   += b5*dN_dy(jj);
              //Klocal(TIp1, TJp1) += b6*dN_dy(jj);

              // convection term

              gradN = grad*(rho*N(jj));

              Db = rho*(vel(0)*dN_dx(jj) + vel(1)*dN_dy(jj));

              gradN(0,0) += Db;
              gradN(1,1) += Db;

              Klocal(TI,   TJ)   += (b8*gradN(0,0));
              Klocal(TI,   TJp1) += (b8*gradN(0,1));
              Klocal(TIp1, TJ)   += (b8*gradN(1,0));
              Klocal(TIp1, TJp1) += (b8*gradN(1,1));

              // pressure term
              Klocal(TI,   TJp2) -= (b1*N(jj));
              Klocal(TIp1, TJp2) -= (b2*N(jj));

              // continuity equation
              Klocal(TIp2, TJ)   += (b8*dN_dx(jj));
              Klocal(TIp2, TJp1) += (b8*dN_dy(jj));

              // SUPG and PSPG stabilisation terms
              fact2 -= ( muTaf*d2N(jj) );

              gradN *= af;

              Dj(0,0) = gradN(0,0) + fact2;
              Dj(0,1) = gradN(0,1);
              Dj(0,2) = dN_dx(jj);
              Dj(1,0) = gradN(1,0);
              Dj(1,1) = gradN(1,1) + fact2;
              Dj(1,2) = dN_dy(jj);

              if(axsy)
              {
                Dj(0,0) -= muTaf*(dN_dx(jj)/rad - N(jj)/rad/rad);
                Dj(1,1) -= muTaf*(dN_dx(jj)/rad);
              }

              // SUPG
              Klocal(TI, TJ)     += Da*Dj(0,0);
              Klocal(TI, TJp1)   += Da*Dj(0,1);
              Klocal(TI, TJp2)   += Da*Dj(0,2);

              Klocal(TIp1, TJ)   += Da*Dj(1,0);
              Klocal(TIp1, TJp1) += Da*Dj(1,1);
              Klocal(TIp1, TJp2) += Da*Dj(1,2);

              Klocal(TI,   TJ)   += ( (tau[0]*af) * b1 * rStab(0) * N(jj) );
              Klocal(TI,   TJp1) += ( (tau[0]*af) * b2 * rStab(0) * N(jj) );
              Klocal(TIp1, TJ)   += ( (tau[0]*af) * b1 * rStab(1) * N(jj) );
              Klocal(TIp1, TJp1) += ( (tau[0]*af) * b2 * rStab(1) * N(jj) );

              // PSPG
              Klocal(TIp2, TJ)   += (b1*Dj(0,0) + b2*Dj(1,0))*tau[1];
              Klocal(TIp2, TJp1) += (b1*Dj(0,1) + b2*Dj(1,1))*tau[1];
              Klocal(TIp2, TJp2) += (b1*Dj(0,2) + b2*Dj(1,2))*tau[1];

              // LSIC stabilisation

              fact2 = rho*af*tau[2];
              Klocal(TI,   TJ)   += (b1*dN_dx(jj))*fact2;
              Klocal(TI,   TJp1) += (b1*dN_dy(jj))*fact2;

              Klocal(TIp1, TJ)   += (b2*dN_dx(jj))*fact2;
              Klocal(TIp1, TJp1) += (b2*dN_dy(jj))*fact2;

              if(axsy)
              {
                  // diffusion term
                  Klocal(TI, TJ)     += (b4 * (mu/rad/rad) * (af*N(jj)) );
                  Klocal(TI, TJp2)   -= (b4 * N(jj)/rad);

                  // continuity equation
                  Klocal(TIp2, TJ)   += (b4 * af*N(jj)/rad);
              }
            }

            Flocal(TI)   -= (b4*res2(0) + b1*stress(0,0) + b2*stress(0,1) );
            Flocal(TIp1) -= (b4*res2(1) + b1*stress(1,0) + b2*stress(1,1) );
            Flocal(TIp2) -= (b4*grad.trace());

            // SUPG stabilisation terms
            Flocal(TI)   -= Da*rStab(0);
            Flocal(TIp1) -= Da*rStab(1);

            // PSPG stabilisation terms
            Flocal(TIp2) -= (tau[1]*(b1*rStab(0)+b2*rStab(1)));

            // LSIC stabilisation terms

            fact2 = tau[2]*rho*grad.trace();

            Flocal(TI)   -= b1*fact2;
            Flocal(TIp1) -= b2*fact2;

            if(axsy)
            {
                Flocal(TI)   -= (b4 * (mu/rad/rad) * vel(0) );
                Flocal(TI)   += (b4 * pres/rad);
                Flocal(TIp2) -= (b4 * vel(0)/rad);
            }
          } // for(ii=0;ii<totnlbf2;ii++)
    }//gp

    return;
}
*/




/*
template<>
void TreeNode<2>::calcStiffnessAndResidualCutFEMFluid(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
    // GFEM for Navier-Stokes
    // with proper stabilisation
    // semi-implicit formulation - type A
    // 
    ///////////////////////////////////////

    int ii, jj, gp, nGauss, tempId, count=0;
    int TI, TIp1, TIp2, TJ, TJp1, TJp2;

    double  JacTemp, Jac, dvol, stabParam, CI=7.0;
    double  fact, fact2, b1, b2, b3, b4, b5, b6, b7, b8;
    double  pres, Da, Db, rad, urdr, urdr2, h2, h, tau[3];

    //double beta[6]; get_stabilisation_beta_wulf(beta);

    VectorXd  NN(totnlbf), dNN_dx(totnlbf), d2NN_dx2(totnlbf), dNN_dy(totnlbf), d2NN_dy2(totnlbf);
    VectorXd  N, dN_dx, d2N_dx2, dN_dy, d2N_dy2, d2N;
    VectorXd  res(3), res2(2), dp(2), Du(2), vel(2), velDot(2), force(2), gradTvel(2), rStab(3);
    VectorXd  velPrev(2), velPrev2(2), velExt(2), velTemp(3);
    MatrixXd  Dj(2, 3), grad(2,2), gradN(2,2), stress(2,2), gradPrev(2,2), matG(3,3);
    myPoint  param, geom;
    Dj.setZero();


    bool   axsy = ((int)elmDat[2] == 1);
    double  rho = elmDat[3];
    double  mu  = elmDat[4];
    double  bforce[2] = {elmDat[5], elmDat[6]};
    double  af = SolnData->td(2);
    double  am = SolnData->td(1);
    double  acceFact = am*SolnData->td(9);
    double  dt = mpapTime.dt;
    double  muTaf = mu*af;

    double *gws;
    myPoint *gps;

    double  hx = bbox.maxBB[0]-bbox.minBB[0];
    double  hy = bbox.maxBB[1]-bbox.minBB[1];

    double volume = hx*hy;

    matG.setZero();
    matG(0,0) = 4.0/hx/hx;
    matG(1,1) = 4.0/hy/hy;

    if(domNums.size() > 1)
    {
      nGauss = Quadrature.gausspoints.size();

      gps = &(Quadrature.gausspoints[0]);
      gws = &(Quadrature.gaussweights[0]);

      JacTemp = 1.0;

      volume = 0.0;
      for(gp=0; gp<nGauss; gp++)
      {
        //cout << gp << '\t' << gws[gp] << endl;
        volume += gws[gp];
      }

      // For 2D problem
      // fact = sqrt(Vc/V); and fact = fact*fact;  ---->  fact = Vc/V;
      // For 3D problem
      // fact = (Vc/V)^(1/3). So, do fact = fact*fact;

      fact = volume/(hx*hy);

      matG(0,0) /= fact;
      matG(1,1) /= fact;
    }
    else
    {
      nGauss = GeomData->gausspoints.size();

      gps = &(GeomData->gausspoints[0]);
      gws = &(GeomData->gaussweights[0]);

      JacTemp = JacMultElem;
    }

    //cout << nGauss << '\t' << volume << endl;
    //printMatrix(matG);
    //h2 = 4.0*volume/PI;
    //h = sqrt(h2);
    //stabParam = h2/(4.0*mu)/degree[0]/degree[1];
    //stabParam = h2/(12.0*mu);
    //stabParam = h2/(12.0*mu);
    //stabParam /= degree[0]/degree[1];
    ////stabParam /= rho;
    //stabParam *= rho;

    //tau[0] = elmDat[8]*stabParam;  // SUPG
    //tau[1] = elmDat[9]*stabParam;  // PSPG
    //tau[2] = elmDat[10]*stabParam; // LSIC
    //cout << tau[0] << '\t' << tau[1] << '\t' << tau[2] << endl;

    count=0;
    for(gp=0; gp<nGauss; gp++)
    {
        param[0]  = 0.5*(knotIncr[0] * gps[gp][0] + knotSum[0]);
        param[1]  = 0.5*(knotIncr[1] * gps[gp][1] + knotSum[1]);

        dvol = gws[gp] * JacTemp;

        geom[0] = GeomData->computeCoord(0, param[0]);
        geom[1] = GeomData->computeCoord(1, param[1]);

          GeomData->computeBasisFunctions2D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy, d2NN_dx2, d2NN_dy2);

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
            //cout << " parent is not NULL " << endl;
            N = SubDivMat*NN;
            dN_dx = SubDivMat*dNN_dx;
            dN_dy = SubDivMat*dNN_dy;
            d2N_dx2 = SubDivMat*d2NN_dx2;
            d2N_dy2 = SubDivMat*d2NN_dy2;
          }

        d2N = d2N_dx2 + d2N_dy2;

        vel.setZero();
        velPrev.setZero();
        velPrev2.setZero();
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

            b1     = SolnData->var1Prev(TI);
            b2     = SolnData->var1Prev(TIp1);

            velPrev(0) += ( b1 * N(ii) );
            velPrev(1) += ( b2 * N(ii) );

            b1     = SolnData->var1Prev2(TI);
            b2     = SolnData->var1Prev2(TIp1);

            velPrev2(0) += ( b1 * N(ii) );
            velPrev2(1) += ( b2 * N(ii) );

            b1     = SolnData->var1Cur(TI);
            b2     = SolnData->var1Cur(TIp1);
            b3     = SolnData->var1(TIp2);

            vel(0) += ( b1 * N(ii) );
            vel(1) += ( b2 * N(ii) );

            grad(0,0) += ( b1 * dN_dx(ii) );
            grad(0,1) += ( b1 * dN_dy(ii) );
            grad(1,0) += ( b2 * dN_dx(ii) );
            grad(1,1) += ( b2 * dN_dy(ii) );

            Du(0)  += ( b1 * d2N(ii) );
            Du(1)  += ( b2 * d2N(ii) );

            pres   += ( b3 * N(ii) );
            dp(0)  += ( b3 * dN_dx(ii) );
            dp(1)  += ( b3 * dN_dy(ii) );

            b1     = SolnData->var1DotCur(TI);
            b2     = SolnData->var1DotCur(TIp1);

            velDot(0) += ( b1 * N(ii) );
            velDot(1) += ( b2 * N(ii) );
          }

          // this is pseudo-stress
          //stress = mu*(grad+grad.transpose());
          stress = mu*grad;
          stress(0,0) -= pres;
          stress(1,1) -= pres;

          force.setZero();

          //force(0) = GeomData->analyDBC->computeForce(0, xx, yy);
          //force(1) = GeomData->analyDBC->computeForce(0, xx, yy);
          //force = 1.0;
          //force = 0.0;
          //cout << force(0) << '\t' << force(1) << endl;

<<<<<<< HEAD:HBsplines/TreeNode9.cpp
          velExt = velPrev;
=======
          //velExt = velPrev;
          velExt = 2.0*velPrev-velPrev2;
>>>>>>> collabchandan:src/HBsplines/TreeNode8.cpp
          //velExt = af*(2.0*velPrev-velPrev2) + (1.0-af)*velPrev;

          gradTvel = grad*velExt;

          res2(0) = rho*(velDot(0) + gradTvel(0) - force(0)) ;
          res2(1) = rho*(velDot(1) + gradTvel(1) - force(1)) ;

          rStab(0) = res2(0) - mu*Du(0) + dp(0) ;
          rStab(1) = res2(1) - mu*Du(1) + dp(1) ;

          if(axsy)
          {
            rad = geom[0];

            urdr  = vel(0)/rad;
            urdr2 = urdr/rad;

            dvol *= (2.0*PI*rad);

            rStab(0) -= mu*(grad(0,0)/rad - urdr2 );
            rStab(1) -= mu*(grad(1,0)/rad );
          }

          // evaluate stabilisation parameters
          //
          velTemp(0) = velPrev(0);
          velTemp(1) = velPrev(1);
          velTemp(2) = 0.0;

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

            b1 = dN_dx[ii]*dvol;
            b2 = dN_dy[ii]*dvol;
            b4 = N[ii]*dvol;

            b5 = muTaf*b1;
            b6 = muTaf*b2;
            b8 = af*b4;

            Da = rho*(velExt(0)*b1 + velExt(1)*b2)*tau[0];

            for(jj=0;jj<totnlbf2;jj++)
            {
              TJ   = ndof*jj;
              TJp1 = TJ+1;
              TJp2 = TJ+2;

              fact2 = rho*acceFact*N(jj);

              // time acceleration term
              fact = b4*fact2 ;

              // diffusion term
              fact += ( b5*dN_dx(jj)+b6*dN_dy(jj) );

              Klocal(TI,   TJ)   += fact;
              Klocal(TIp1, TJp1) += fact;

              // convection term - semi-implicit type A

              Db = rho*(velExt(0)*dN_dx(jj) + velExt(1)*dN_dy(jj));

              Klocal(TI,   TJ)   += (b8*Db);
              Klocal(TIp1, TJp1) += (b8*Db);

              // pressure term
              Klocal(TI,   TJp2) -= (b1*N(jj));
              Klocal(TIp1, TJp2) -= (b2*N(jj));

              // continuity equation
              Klocal(TIp2, TJ)   += (b8*dN_dx(jj));
              Klocal(TIp2, TJp1) += (b8*dN_dy(jj));

              // SUPG and PSPG stabilisation terms
              fact2 -= ( muTaf*d2N(jj) );

              Db *= af;

              Dj(0,0) = Db + fact2;
              Dj(0,1) = 0.0;
              Dj(0,2) = dN_dx(jj);
              Dj(1,0) = 0.0;
              Dj(1,1) = Db + fact2;
              Dj(1,2) = dN_dy(jj);

              if(axsy)
              {
                Dj(0,0) -= muTaf*(dN_dx(jj)/rad - N(jj)/rad/rad);
                Dj(1,1) -= muTaf*(dN_dx(jj)/rad);
              }

              // SUPG
              Klocal(TI, TJ)     += Da*Dj(0,0);
              Klocal(TI, TJp1)   += Da*Dj(0,1);
              Klocal(TI, TJp2)   += Da*Dj(0,2);

              Klocal(TIp1, TJ)   += Da*Dj(1,0);
              Klocal(TIp1, TJp1) += Da*Dj(1,1);
              Klocal(TIp1, TJp2) += Da*Dj(1,2);

              // PSPG
              Klocal(TIp2, TJ)   += (b1*Dj(0,0) + b2*Dj(1,0))*tau[1];
              Klocal(TIp2, TJp1) += (b1*Dj(0,1) + b2*Dj(1,1))*tau[1];
              Klocal(TIp2, TJp2) += (b1*Dj(0,2) + b2*Dj(1,2))*tau[1];

              // LSIC stabilisation

              fact2 = rho*tau[2];
              Klocal(TI,   TJ)   += (b1*dN_dx(jj))*fact2;
              Klocal(TI,   TJp1) += (b1*dN_dy(jj))*fact2;

              Klocal(TIp1, TJ)   += (b2*dN_dx(jj))*fact2;
              Klocal(TIp1, TJp1) += (b2*dN_dy(jj))*fact2;

              if(axsy)
              {
                  // diffusion term
                  Klocal(TI, TJ)     += (b4 * (mu/rad/rad) * (af*N(jj)) );
                  Klocal(TI, TJp2)   -= (b4 * N(jj)/rad);

                  // continuity equation
                  Klocal(TIp2, TJ)   += (b4 * af*N(jj)/rad);
              }
            }

            Flocal(TI)   -= (b4*res2(0) + b1*stress(0,0) + b2*stress(0,1) );
            Flocal(TIp1) -= (b4*res2(1) + b1*stress(1,0) + b2*stress(1,1) );
            Flocal(TIp2) -= (b4*grad.trace());

            // SUPG stabilisation terms
            Flocal(TI)   -= Da*rStab(0);
            Flocal(TIp1) -= Da*rStab(1);

            // PSPG stabilisation terms
            Flocal(TIp2) -= (tau[1]*(b1*rStab(0)+b2*rStab(1)));

            // LSIC stabilisation terms

            fact2 = rho*tau[2]*grad.trace();
            Flocal(TI)   -= b1*fact2;
            Flocal(TIp1) -= b2*fact2;

            if(axsy)
            {
                Flocal(TI)   -= (b4 * (mu/rad/rad) * vel(0) );
                Flocal(TI)   += (b4 * pres/rad);
                Flocal(TIp2) -= (b4 * vel(0)/rad);
            }
          } // for(ii=0;ii<totnlbf2;ii++)
      //} //if(within)
    }//gp

    return;
}
*/



//
template<>
void TreeNode<2>::calcStiffnessAndResidualCutFEMFluid(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
    // GFEM for Navier-Stokes
    // with proper stabilisation
    // SUPG/PSPG stabilisation
    // semi-implicit formulation - type B
    // 
    ///////////////////////////////////////

    int  ii=0, jj=0, gp=0, nGauss=0, tempId=0, count=0;
    int  TI=0, TIp1=0, TIp2=0, TJ=0, TJp1=0, TJp2=0;

    double  JacTemp=0.0, Jac=0.0, dvol=0.0, stabParam=0.0, CI=7.0;
    double  fact=0.0, fact2=0.0, b1=0.0, b2=0.0, b3=0.0, b4=0.0, b5=0.0, b6=0.0, b7=0.0, b8=0.0;
    double  pres=0.0, Da=0.0, Db=0.0, rad=0.0, urdr=0.0, urdr2=0.0, h2=0.0, h=0.0, tau[3];

    //double beta[6]; get_stabilisation_beta_wulf(beta);

    VectorXd  NN(totnlbf), dNN_dx(totnlbf), d2NN_dx2(totnlbf), dNN_dy(totnlbf), d2NN_dy2(totnlbf);
    VectorXd  N(totnlbf), dN_dx(totnlbf), d2N_dx2(totnlbf), dN_dy(totnlbf), d2N_dy2(totnlbf), d2N(totnlbf);
    VectorXd  res(3), res2(2), dp(2), Du(2), vel(2), velDot(2), velTemp(3);
    VectorXd  force(2), gradTvel(2), rStab(3), velPrev(2);
    MatrixXd  Dj(2, 3), grad(2,2), gradN(2,2), stress(2,2), gradPrev(2,2), matG(3,3);
    myPoint  param, geom;
    Dj.setZero();


    bool   axsy = ((int)elmDat[2] == 1);
    double  rho = elmDat[3];
    double  mu  = elmDat[4];
    double  bforce[2] = {elmDat[5], elmDat[6]};
    double  af = SolnData->td(2);
    double  am = SolnData->td(1);
    double  acceFact = am*SolnData->td(9);
    double  dt = mpapTime.dt;
    double  muTaf = mu*af;

    double *gws;
    myPoint *gps;

    double  hx = bbox.maxBB[0]-bbox.minBB[0];
    double  hy = bbox.maxBB[1]-bbox.minBB[1];

    double volume = hx*hy;

    matG.setZero();
    matG(0,0) = 4.0/hx/hx;
    matG(1,1) = 4.0/hy/hy;

    if(domNums.size() > 1)
    {
      nGauss = Quadrature.gausspoints.size();

      gps = &(Quadrature.gausspoints[0]);
      gws = &(Quadrature.gaussweights[0]);

      JacTemp = 1.0;

      volume = 0.0;
      for(gp=0; gp<nGauss; gp++)
        volume += gws[gp];

      // For 2D problem
      // fact = sqrt(Vc/V); and fact = fact*fact;  ---->  fact = Vc/V;
      // For 3D problem
      // fact = (Vc/V)^(1/3). So, do fact = fact*fact;

      fact = volume/(hx*hy);

      matG(0,0) /= fact;
      matG(1,1) /= fact;

      //matG(0,0) = 4.0/0.909/0.909;
      //matG(1,1) = 4.0/0.380952/0.380952;
    }
    else
    {
      nGauss = GeomData->gausspoints.size();

      gps = &(GeomData->gausspoints[0]);
      gws = &(GeomData->gaussweights[0]);

      JacTemp = JacMultElem;
    }

    //cout << volume << '\t' << volume << endl;
    //h2 = 4.0*volume/PI;
    //stabParam = h2/(12.0*mu);
    ////stabParam /= rho;
    //stabParam *= rho;

    //tau[0] = elmDat[8]*stabParam;  // SUPG
    //tau[1] = elmDat[9]*stabParam;  // PSPG
    //tau[2] = elmDat[10]*stabParam; // LSIC
    //cout << tau[0] << '\t' << tau[1] << '\t' << tau[2] << endl;

    count=0;
    for(gp=0; gp<nGauss; gp++)
    {
        param[0]  = 0.5*(knotIncr[0] * gps[gp][0] + knotSum[0]);
        param[1]  = 0.5*(knotIncr[1] * gps[gp][1] + knotSum[1]);

        dvol = gws[gp] * JacTemp;

        geom[0] = GeomData->computeCoord(0, param[0]);
        geom[1] = GeomData->computeCoord(1, param[1]);

        if(parent == NULL)  //cout << " parent is NULL " << endl;
        {
          GeomData->computeBasisFunctions2D(knotBegin, knotIncr, param, N, dN_dx, dN_dy, d2N_dx2, d2N_dy2);
        }
        else
        {
          GeomData->computeBasisFunctions2D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy, d2NN_dx2, d2NN_dy2);
          //cout << " parent is not NULL " << endl;
          N = SubDivMat*NN;
          dN_dx = SubDivMat*dNN_dx;
          dN_dy = SubDivMat*dNN_dy;
          d2N_dx2 = SubDivMat*d2NN_dx2;
          d2N_dy2 = SubDivMat*d2NN_dy2;
        }

        d2N = d2N_dx2 + d2N_dy2;

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

            b1     = SolnData->var1Prev(TI);
            b2     = SolnData->var1Prev(TIp1);

            velPrev(0) += ( b1 * N(ii) );
            velPrev(1) += ( b2 * N(ii) );

            gradPrev(0,0) += ( b1 * dN_dx(ii) );
            gradPrev(0,1) += ( b1 * dN_dy(ii) );
            gradPrev(1,0) += ( b2 * dN_dx(ii) );
            gradPrev(1,1) += ( b2 * dN_dy(ii) );

            b1     = SolnData->var1Cur(TI);
            b2     = SolnData->var1Cur(TIp1);
            b3     = SolnData->var1(TIp2);

            vel(0) += ( b1 * N(ii) );
            vel(1) += ( b2 * N(ii) );

            grad(0,0) += ( b1 * dN_dx(ii) );
            grad(0,1) += ( b1 * dN_dy(ii) );
            grad(1,0) += ( b2 * dN_dx(ii) );
            grad(1,1) += ( b2 * dN_dy(ii) );

            Du(0)  += ( b1 * d2N(ii) );
            Du(1)  += ( b2 * d2N(ii) );

            pres   += ( b3 * N(ii) );
            dp(0)  += ( b3 * dN_dx(ii) );
            dp(1)  += ( b3 * dN_dy(ii) );

            b1     = SolnData->var1DotCur(TI);
            b2     = SolnData->var1DotCur(TIp1);

            velDot(0) += ( b1 * N(ii) );
            velDot(1) += ( b2 * N(ii) );
          }

          //stress = mu*(grad+grad.transpose());
          // this is pseudo-stress
          stress = mu*grad;
          stress(0,0) -= pres;
          stress(1,1) -= pres;

          force.setZero();

          //force(0) = GeomData->analyDBC->computeForce(0, xx, yy);
          //force(1) = GeomData->analyDBC->computeForce(0, xx, yy);
          //force = 1.0;
          //force = 0.0;
          //cout << force(0) << '\t' << force(1) << endl;

          //force[0] -= bforce[0]*timeFunction[0].prop;
          //force[1] -= bforce[1]*timeFunction[0].prop;

          force[0] -= bforce[0];
          force[1] -= bforce[1];

          //gradTvel = grad*velPrev;
          gradTvel = gradPrev*vel + grad*velPrev - gradPrev*velPrev;

          res2(0) = rho*(velDot(0) + gradTvel(0) - force(0)) ;
          res2(1) = rho*(velDot(1) + gradTvel(1) - force(1)) ;

          rStab(0) = res2(0) - mu*Du(0) + dp(0) ;
          rStab(1) = res2(1) - mu*Du(1) + dp(1) ;

          if(axsy)
          {
            rad = geom[0];

            urdr  = vel(0)/rad;
            urdr2 = urdr/rad;

            dvol *= (2.0*PI*rad);

            rStab(0) -= mu*(grad(0,0)/rad - urdr2 );
            rStab(1) -= mu*(grad(1,0)/rad );
          }

          // evaluate stabilisation parameters
          //
          velTemp(0) = velPrev(0);
          velTemp(1) = velPrev(1);
          velTemp(2) = 0.0;

          //evaluateStabParams_algo1(&velTemp(0), h, rho, mu, dt,  beta, tau);
          //evaluateStabParams_algo2(&velTemp(0), h, rho, mu, dt,  beta, tau);
          evaluateStabParams_algo3(velTemp, matG, dt, rho, mu, CI, tau);

          tau[0] *= elmDat[8];  // SUPG
          tau[1] *= elmDat[9];  // PSPG
          tau[2] *= elmDat[10]; // LSIC

          //cout << " tau = " << tau[0] << '\t' << tau[1] << '\t' << tau[2] << endl;

          for(ii=0;ii<totnlbf2;ii++)
          {
            TI   = ndof*ii;
            TIp1 = TI+1;
            TIp2 = TI+2;

            b1 = dN_dx[ii]*dvol;
            b2 = dN_dy[ii]*dvol;
            b4 = N[ii]*dvol;

            b5 = muTaf*b1;
            b6 = muTaf*b2;
            b8 = af*b4;

            Da = rho*(velPrev(0)*b1 + velPrev(1)*b2)*tau[0];

            for(jj=0;jj<totnlbf2;jj++)
            {
              TJ   = ndof*jj;
              TJp1 = TJ+1;
              TJp2 = TJ+2;

              fact2 = rho*acceFact*N(jj);

              // time acceleration term
              fact = b4*fact2 ;

              // diffusion term
              fact += ( b5*dN_dx(jj)+b6*dN_dy(jj) );

              Klocal(TI,   TJ)   += fact;
              Klocal(TIp1, TJp1) += fact;

              // convection term - semi-implicit type B

              gradN = gradPrev*(rho*N(jj));
              //gradN.setZero();

              Db = rho*(velPrev(0)*dN_dx(jj) + velPrev(1)*dN_dy(jj));

              gradN(0,0) += Db;
              gradN(1,1) += Db;

              Klocal(TI,   TJ)   += (b8*gradN(0,0));
              Klocal(TI,   TJp1) += (b8*gradN(0,1));
              Klocal(TIp1, TJ)   += (b8*gradN(1,0));
              Klocal(TIp1, TJp1) += (b8*gradN(1,1));

              // pressure term
              Klocal(TI,   TJp2) -= (b1*N(jj));
              Klocal(TIp1, TJp2) -= (b2*N(jj));

              // continuity equation
              Klocal(TIp2, TJ)   += (b8*dN_dx(jj));
              Klocal(TIp2, TJp1) += (b8*dN_dy(jj));

              // SUPG and PSPG stabilisation terms
              fact2 -= ( muTaf*d2N(jj) );

              gradN *= af;
              //gradN.setZero();

              Dj(0,0) = gradN(0,0) + fact2;
              Dj(0,1) = gradN(0,1);
              Dj(0,2) = dN_dx(jj);
              Dj(1,0) = gradN(1,0);
              Dj(1,1) = gradN(1,1) + fact2;
              Dj(1,2) = dN_dy(jj);

              if(axsy)
              {
                Dj(0,0) -= muTaf*(dN_dx(jj)/rad - N(jj)/rad/rad);
                Dj(1,1) -= muTaf*(dN_dx(jj)/rad);
              }

              // SUPG
              Klocal(TI, TJ)     += Da*Dj(0,0);
              Klocal(TI, TJp1)   += Da*Dj(0,1);
              Klocal(TI, TJp2)   += Da*Dj(0,2);

              Klocal(TIp1, TJ)   += Da*Dj(1,0);
              Klocal(TIp1, TJp1) += Da*Dj(1,1);
              Klocal(TIp1, TJp2) += Da*Dj(1,2);

              // PSPG
              Klocal(TIp2, TJ)   += (b1*Dj(0,0) + b2*Dj(1,0))*tau[1];
              Klocal(TIp2, TJp1) += (b1*Dj(0,1) + b2*Dj(1,1))*tau[1];
              Klocal(TIp2, TJp2) += (b1*Dj(0,2) + b2*Dj(1,2))*tau[1];

              // LSIC stabilisation

              fact2 = rho*af*tau[2];
              //fact2 = af*tau[2];

              Klocal(TI,   TJ)   += (b1*dN_dx(jj))*fact2;
              Klocal(TI,   TJp1) += (b1*dN_dy(jj))*fact2;

              Klocal(TIp1, TJ)   += (b2*dN_dx(jj))*fact2;
              Klocal(TIp1, TJp1) += (b2*dN_dy(jj))*fact2;

              if(axsy)
              {
                  // diffusion term
                  Klocal(TI, TJ)     += (b4 * (mu/rad/rad) * (af*N(jj)) );
                  Klocal(TI, TJp2)   -= (b4 * N(jj)/rad);

                  // continuity equation
                  Klocal(TIp2, TJ)   += (b4 * af*N(jj)/rad);
              }
            }

            Flocal(TI)   -= (b4*res2(0) + b1*stress(0,0) + b2*stress(0,1) );
            Flocal(TIp1) -= (b4*res2(1) + b1*stress(1,0) + b2*stress(1,1) );
            Flocal(TIp2) -= (b4*grad.trace());

            // SUPG stabilisation terms
            Flocal(TI)   -= Da*rStab(0);
            Flocal(TIp1) -= Da*rStab(1);

            // PSPG stabilisation terms
            Flocal(TIp2) -= (tau[1]*(b1*rStab(0)+b2*rStab(1)));

            // LSIC stabilisation terms
            fact2 = tau[2]*rho*grad.trace();
            //fact2 = tau[2]*grad.trace();

            Flocal(TI)   -= b1*fact2;
            Flocal(TIp1) -= b2*fact2;

            if(axsy)
            {
                Flocal(TI)   -= (b4 * (mu/rad/rad) * vel(0) );
                Flocal(TI)   += (b4 * pres/rad);
                Flocal(TIp2) -= (b4 * vel(0)/rad);
            }
          } // for(ii=0;ii<totnlbf2;ii++)
      //} //if(within)
    }//gp

    return;
}
//




/*
template<>
void TreeNode<2>::calcStiffnessAndResidualCutFEMFluid(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
    // GFEM for Navier-Stokes
    // with proper stabilisation
    // GLS stabilisation
    // semi-implicit formulation - type B
    // 
    ///////////////////////////////////////

    int ii, jj, gp, nGauss, tempId, count=0;
    int TI, TIp1, TIp2, TJ, TJp1, TJp2;

    double  JacTemp, Jac, dvol, stabParam;
    double  fact, fact2, b1, b2, b3, b4, b5, b6, b7, b8, acceFact, dt;
    double  pres, Da, Db, af, am, muTaf, rad, urdr, urdr2, h2, h, tau[3];

    //double beta[6]; get_stabilisation_beta_wulf(beta);

    VectorXd  NN(totnlbf), dNN_dx(totnlbf), d2NN_dx2(totnlbf), dNN_dy(totnlbf), d2NN_dy2(totnlbf);
    VectorXd  N, dN_dx, d2N_dx2, dN_dy, d2N_dy2, d2N;
    VectorXd  res(3), res2(2), dp(2), Du(2), vel(2), velDot(2);
    VectorXd  force(2), gradTvel(2), rStab(3), velPrev(2);
    MatrixXd  D(forAssyVec.size(), 3), grad(2,2), gradN(2,2), stress(2,2), gradPrev(2,2);
    myPoint  param, geom, velTemp;
    D.setZero();

    bool   axsy = ((int)elmDat[2] == 1);
    double  rho = elmDat[3];
    double  mu  = elmDat[4];

    af = SolnData->td(2);
    am = SolnData->td(1);
    acceFact = am*SolnData->td(9);
    dt = mpapTime.dt;

    muTaf = mu*af;

    double *gws;
    myPoint *gps;
    
    if(domNums.size() > 1)
    {
      nGauss = Quadrature.gausspoints.size();
      
      gps = &(Quadrature.gausspoints[0]);
      gws = &(Quadrature.gaussweights[0]);
      
      JacTemp = 1.0;

      volume = 0.0;
      for(gp=0; gp<nGauss; gp++)
        volume += gws[gp];
    }
    else
    {
      nGauss = GeomData->gausspoints.size();

      gps = &(GeomData->gausspoints[0]);
      gws = &(GeomData->gaussweights[0]);
      
      JacTemp = JacMultElem;
    }

    //cout << " nGauss " << nGauss << '\t' << tempId << endl;

//    cout << volume << '\t' << volume << endl;
    h2 = 4.0*volume/PI;
    h = sqrt(h2);

    stabParam = h2/(12.0*mu);
    stabParam /= degree[0]/degree[1];
    //stabParam /= rho;
    stabParam *= rho;

    tau[0] = elmDat[8]*stabParam;//rho;      // SUPG
    tau[1] = elmDat[9]*stabParam;  // PSPG
    tau[2] = elmDat[10]*stabParam; // LSIC

    //cout << tau[0] << '\t' << tau[1] << '\t' << tau[2] << endl;

    count=0;
    for(gp=0; gp<nGauss; gp++)
    {
        param[0]  = 0.5*(knotIncr[0] * gps[gp][0] + knotSum[0]);
        param[1]  = 0.5*(knotIncr[1] * gps[gp][1] + knotSum[1]);

        dvol = gws[gp] * JacTemp;

        geom[0] = GeomData->computeCoord(0, param[0]);
        geom[1] = GeomData->computeCoord(1, param[1]);

        //cout << uu << '\t' << vv << endl;
        //cout << xx << '\t' << yy << endl;

          GeomData->computeBasisFunctions2D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy, d2NN_dx2, d2NN_dy2);

          if(parent == NULL)
          {
            //cout << " parent is NULL " << endl;
            N = NN;
            dN_dx = dNN_dx;
            dN_dy = dNN_dy;
            d2N_dx2 = d2NN_dx2;
            d2N_dy2 = d2NN_dy2;
          }
          else
          {
            //cout << " parent is not NULL " << endl;
            N = SubDivMat*NN;
            dN_dx = SubDivMat*dNN_dx;
            dN_dy = SubDivMat*dNN_dy;
            d2N_dx2 = SubDivMat*d2NN_dx2;
            d2N_dy2 = SubDivMat*d2NN_dy2;
          }

        //cout <<  " gpDomainId = " << gpDomainId << '\t' << totnlbf2 << endl;
        //printVector(GlobalBasisFuncs);

        d2N = d2N_dx2 + d2N_dy2;

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

            b1     = SolnData->var1Prev(TI);
            b2     = SolnData->var1Prev(TIp1);

            velPrev(0) += ( b1 * N(ii) );
            velPrev(1) += ( b2 * N(ii) );

            gradPrev(0,0) += ( b1 * dN_dx(ii) );
            gradPrev(0,1) += ( b1 * dN_dy(ii) );
            gradPrev(1,0) += ( b2 * dN_dx(ii) );
            gradPrev(1,1) += ( b2 * dN_dy(ii) );

            b1     = SolnData->var1Cur(TI);
            b2     = SolnData->var1Cur(TIp1);
            b3     = SolnData->var1(TIp2);

            vel(0) += ( b1 * N(ii) );
            vel(1) += ( b2 * N(ii) );

            grad(0,0) += ( b1 * dN_dx(ii) );
            grad(0,1) += ( b1 * dN_dy(ii) );
            grad(1,0) += ( b2 * dN_dx(ii) );
            grad(1,1) += ( b2 * dN_dy(ii) );

            Du(0)  += ( b1 * d2N(ii) );
            Du(1)  += ( b2 * d2N(ii) );

            pres   += ( b3 * N(ii) );
            dp(0)  += ( b3 * dN_dx(ii) );
            dp(1)  += ( b3 * dN_dy(ii) );

            b1     = SolnData->var1DotCur(TI);
            b2     = SolnData->var1DotCur(TIp1);

            velDot(0) += ( b1 * N(ii) );
            velDot(1) += ( b2 * N(ii) );
          }

          // this is pseudo-stress
          //stress = mu*(grad+grad.transpose());
          stress = mu*grad;
          stress(0,0) -= pres;
          stress(1,1) -= pres;

          force.setZero();

          //force(0) = GeomData->analyDBC->computeForce(0, xx, yy);
          //force(1) = GeomData->analyDBC->computeForce(0, xx, yy);
          //force = 1.0;
          //force = 0.0;
          //cout << force(0) << '\t' << force(1) << endl;

          gradTvel = gradPrev*vel + grad*velPrev - gradPrev*velPrev;

          res2(0) = rho*(velDot(0) + gradTvel(0) - force(0)) ;
          res2(1) = rho*(velDot(1) + gradTvel(1) - force(1)) ;

          rStab(0) = res2(0) - mu*Du(0) + dp(0) ;
          rStab(1) = res2(1) - mu*Du(1) + dp(1) ;

          if(axsy)
          {
            rad = geom[0];

            urdr  = vel(0)/rad;
            urdr2 = urdr/rad;

            dvol *= (2.0*PI*rad);
            
            rStab(0) -= mu*(grad(0,0)/rad - urdr2 );
            rStab(1) -= mu*(grad(1,0)/rad );
          }

          // evaluate stabilisation parameters
          //
          velTemp(0) = velPrev(0);
          velTemp(1) = velPrev(1);
          velTemp(2) = 0.0;
          eval_tau(&velTemp(0), h, rho, mu, dt,  beta, tau);

          tau[0] = elmDat[8]*stabParam;  // SUPG
          tau[1] = elmDat[9]*stabParam;  // PSPG
          tau[2] = elmDat[10]*stabParam; // LSIC

          //tau[0] *= rho;
          //tau[1] /= rho;

          //cout << " tau = " << tau[0] << '\t' << tau[1] << '\t' << tau[2] << endl;

          for(ii=0;ii<totnlbf2;ii++)
          {
            TI   = ndof*ii;
            TIp1 = TI+1;
            TIp2 = TI+2;

            b1 = dN_dx[ii]*dvol;
            b2 = dN_dy[ii]*dvol;
            b4 = N[ii]*dvol;

            b5 = muTaf*b1;
            b6 = muTaf*b2;
            b8 = af*b4;

            // GLS stabilisation terms

            fact = rho*acceFact*N[ii] + af*(rho*(velPrev(0)*dN_dx[ii] + velPrev(1)*dN_dy[ii]) - mu*d2N[ii]);

            D(TI,0)   = fact ;
            D(TIp1,0) = 0.0;
            D(TIp2,0) = dN_dx[ii];

            D(TI,1)   = 0.0;
            D(TIp1,1) = fact;
            D(TIp2,1) = dN_dy[ii];

            D(TI,2)   = af*dN_dx[ii];
            D(TIp1,2) = af*dN_dy[ii];
            D(TIp2,2) = 0.0;


            for(jj=0;jj<totnlbf2;jj++)
            {
              TJ   = ndof*jj;
              TJp1 = TJ+1;
              TJp2 = TJ+2;

              fact2 = rho*acceFact*N(jj);

              // time acceleration term
              fact = b4*fact2 ;

              // diffusion term
              fact += ( b5*dN_dx(jj)+b6*dN_dy(jj) );

              Klocal(TI,   TJ)   += fact;
              Klocal(TIp1, TJp1) += fact;

              // convection term - semi-implicit type B

              gradN = gradPrev*(rho*N(jj));

              Db = rho*(velPrev(0)*dN_dx(jj) + velPrev(1)*dN_dy(jj));

              gradN(0,0) += Db;
              gradN(1,1) += Db;

              Klocal(TI,   TJ)   += (b8*gradN(0,0));
              Klocal(TI,   TJp1) += (b8*gradN(0,1));
              Klocal(TIp1, TJ)   += (b8*gradN(1,0));
              Klocal(TIp1, TJp1) += (b8*gradN(1,1));

              // pressure term
              Klocal(TI,   TJp2) -= (b1*N(jj));
              Klocal(TIp1, TJp2) -= (b2*N(jj));

              // continuity equation
              Klocal(TIp2, TJ)   += (b8*dN_dx(jj));
              Klocal(TIp2, TJp1) += (b8*dN_dy(jj));

              if(axsy)
              {
                  // diffusion term
                  Klocal(TI, TJ)     += (b4 * (mu/rad/rad) * (af*N(jj)) );
                  Klocal(TI, TJp2)   -= (b4 * N(jj)/rad);

                  // continuity equation
                  Klocal(TIp2, TJ)   += (b4 * af*N(jj)/rad);
              }
            }

            Flocal(TI)   -= (b4*res2(0) + b1*stress(0,0) + b2*stress(0,1) );
            Flocal(TIp1) -= (b4*res2(1) + b1*stress(1,0) + b2*stress(1,1) );
            Flocal(TIp2) -= (b4*grad.trace());

            if(axsy)
            {
                Flocal(TI)   -= (b4 * (mu/rad/rad) * vel(0) );
                Flocal(TI)   += (b4 * pres/rad);
                Flocal(TIp2) -= (b4 * vel(0)/rad);
            }
          } // for(ii=0;ii<totnlbf2;ii++)

          // GLS stabilisation terms

          res.setZero();
          //res(0) = dt*computeForceCur(0, N) ;
          //res(1) = dt*computeForceCur(1, N) ;

          res(0) -= rStab(0) ;
          res(1) -= rStab(1) ;
          res(2) -= grad.trace() ;

          dvol *= tau[1];

          Klocal += ((dvol*D)*D.transpose());

          Flocal += (D*(dvol*res));

      //} //if(within)
    }//gp

    return;
}
*/



template<>
void TreeNode<2>::applyBoundaryConditionsAtApointCutFEMFluid(myDataIntegrateCutFEM& myData)
{
  // compute stiffness and force vectors corresponding 
  // to Nitsche method of applying interface conditions
  // 
  // diagonal terms
  
    int ii, jj, TI, TIp1, TIp2, TJ, TJp1, TJp2;

    double  fact, fact1, fact2, res, trac;
    double  Ta1, Ta2, Tb1, Tb2, mu1, mu2, specVal;
    double  bb1, bb2, u1, u2, t1, t2, pres;

    VectorXd  NN(totnlbf), dNN_dx(totnlbf), dNN_dy(totnlbf), N, dN_dx, dN_dy;
    VectorXd  vel1(2), vel2(2), trac1(2), trac2(2);
    MatrixXd  stress1(2,2), stress2(2,2);
    myPoint   normal1, normal2;

    bool   axsy = ((int)elmDat[2] == 1);
    double  rho = elmDat[3];
    double  af = SolnData->td(2);

    GeomData->computeBasisFunctions2D(knotBegin, knotIncr, myData.param, NN, dNN_dx, dNN_dy);

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

    // at this moment no jumps in velocity or tractions is considered

    vel1[0] = myData.specVal[0];
    vel1[1] = myData.specVal[1];
    vel2 = vel1;

    vel1(0) -= computeValueCur(0, N);
    vel1(1) -= computeValueCur(1, N);

    vel2(0) -= computeValue2Cur(0, N);
    vel2(1) -= computeValue2Cur(1, N);

    ///////////////////////////////////
    // compute the stresses and tractions
    ///////////////////////////////////

    // domain 1

    mu1 = elmDat[4];
    normal1 = -myData.normal;

    stress1(0,0) = mu1*computeValueCur(0, dN_dx);
    stress1(0,1) = mu1*computeValueCur(0, dN_dy);
    stress1(1,0) = mu1*computeValueCur(1, dN_dx);
    stress1(1,1) = mu1*computeValueCur(1, dN_dy);
    pres         = computeValue(2, N);

    stress1(0,0) -= pres;
    stress1(1,1) -= pres;

    trac1[0] = 0.0 + ( stress1(0,0)*normal1[0] + stress1(0,1)*normal1[1] );
    trac1[1] = 0.0 + ( stress1(1,0)*normal1[0] + stress1(1,1)*normal1[1] );

    // domain 2

    mu2 = elmDat[5];
    normal2 = -normal1;

    stress2(0,0) = mu2*computeValue2Cur(0, dN_dx);
    stress2(0,1) = mu2*computeValue2Cur(0, dN_dy);
    stress2(1,0) = mu2*computeValue2Cur(1, dN_dx);
    stress2(1,1) = mu2*computeValue2Cur(1, dN_dy);
    pres         = computeValue2(2, N);

    stress2(0,0) -= pres;
    stress2(1,1) -= pres;

    trac2[0] = 0.0 + ( stress2(0,0)*normal2[0] + stress1(0,1)*normal2[1] );
    trac2[1] = 0.0 + ( stress2(1,0)*normal2[0] + stress1(1,1)*normal2[1] );

    //for(ii=0;ii<totnlbf;ii++)
      //printf(" \t %14.8f \n", N[ii]);

    //specVal = GeomData->analyDBC->computeValue(0, myData.geom[0], myData.geom[1]);

    //cout << " myData.PENALTY = " << myData.PENALTY << endl;

    for(ii=0;ii<totnlbf2;ii++)
    {
        TI   = ndof*ii;
        TIp1 = TI+1;
        TIp2 = TI+2;

        bb1 = N[ii] * myData.dvol;
        bb2 = bb1 * myData.PENALTY ;

        Ta1 = (myData.dvol*mu1)*( dN_dx(ii)*normal1[0] + dN_dy(ii)*normal1[1] );

        for(jj=0;jj<totnlbf2;jj++)
        {
          TJ   = ndof*jj;
          TJp1 = TJ+1;
          TJp2 = TJ+2;

          fact = af* bb2 * N[jj];
          // stabilisation term
          myData.K1(TI,   TJ)   += fact;
          myData.K1(TIp1, TJp1) += fact;

          // Nitsche terms
          Tb1 = af*mu1*( dN_dx(jj)*normal1[0] + dN_dy(jj)*normal1[1] );

          fact1 = -normal1[0]*N(jj);
          fact2 = -normal1[1]*N(jj);

          myData.K1(TI, TJ)      -= (bb1*Tb1);
          myData.K1(TI, TJp2)    -= (bb1*fact1);

          myData.K1(TI, TJ)      -= (Ta1*af*N(jj))*myData.NitscheFact;
          myData.K1(TIp2, TJ)    -= (bb1*fact1*af*myData.NitscheFact);

          myData.K1(TIp1, TJp1)  -= (bb1*Tb1);
          myData.K1(TIp1, TJp2)  -= (bb1*fact2);

          myData.K1(TIp1, TJp1)  -= (Ta1*af*N(jj))*myData.NitscheFact;
          myData.K1(TIp2, TJp1)  -= (bb1*fact2*af*myData.NitscheFact);
        }

        // stabilisation terms
        myData.F1(TI)   += (bb2*vel1[0]);
        myData.F1(TIp1) += (bb2*vel1[1]);

        // Nitsche terms

        myData.F1(TI)   -= (-bb1*trac1[0]);
        myData.F1(TI)   -= (Ta1*vel1[0])*myData.NitscheFact;
        myData.F1(TIp2) -= (bb1*(-normal1[0]*vel1[0]))*myData.NitscheFact;

        myData.F1(TIp1) -= (-bb1*trac1[1]);
        myData.F1(TIp1) -= (Ta1*vel1[1])*myData.NitscheFact;
        myData.F1(TIp2) -= (bb1*(-normal1[1]*vel1[1]))*myData.NitscheFact;
        
        // terms for the monolithic scheme

        //myData.K2(0, TIp1)  -= (af*bb2);
        myData.K2(0, TIp1)  += (af*Ta1);
        myData.K2(0, TIp2)  -= (bb1*normal1[1]);

        //myData.Kc(0, 0)     += (af*myData.dvol*myData.PENALTY);

        //myData.F2(0)   -= (myData.dvol*myData.PENALTY*vel1[1]);
        myData.F2(0)   -= (myData.dvol*trac1[1]);
    } // for(ii=0;ii<totnlbf2;ii++)
    //

   return;
}



template<>
void TreeNode<2>::applyGhostPenaltyCutFEM(myDataIntegrateCutFEM& myData)
{

   return;
}



template<>
void TreeNode<2>::applyBoundaryConditionsAtApointCutFEMFluid2(myDataIntegrateCutFEM& myData)
{
  // compute stiffness and force vectors corresponding 
  // to Nitsche method of applying interface conditions
  // 
  // coupling terms


   return;
}



template<>
int TreeNode<2>::computeGaussPointsAdapIntegrationBoundary(int side, int refLev1, int refLev2, int inclFlag, int domainCur)
{
    myPoint  pt1, pt2, ptTemp;

    //vector<myPoint>  myLine;

    int  domTemp, ii, jj, kk, nGauss;

    double  val, param, dist;
    bool  flag;

    switch(side)
    {
        case 0:

          pt1[0] = bbox.minBB[0];   pt1[1] = bbox.minBB[1];
          pt2[0] = bbox.minBB[0];   pt2[1] = bbox.maxBB[1];
          param  = -1.0;

        break;

        case 1:

          pt1[0] = bbox.maxBB[0];   pt1[1] = bbox.minBB[1];
          pt2[0] = bbox.maxBB[0];   pt2[1] = bbox.maxBB[1];
          param  = 1.0;

        break;

        case 2:

          pt1[0] = bbox.minBB[0];   pt1[1] = bbox.minBB[1];
          pt2[0] = bbox.maxBB[0];   pt2[1] = bbox.minBB[1];
          param  = -1.0;

        break;

        case 3:

          pt1[0] = bbox.minBB[0];   pt1[1] = bbox.maxBB[1];
          pt2[0] = bbox.maxBB[0];   pt2[1] = bbox.maxBB[1];
          param  = 1.0;

        break;

        default :

          cout << " Invalid 'side' value in GeomDataHBSplines::getBoundaryGPs2D " << endl;
        break;
    } //switch(side)

    // identify the edge of the subtriangle that is on the boundary of the background grid

  /*
  flag = false;
  for(vector<myPoly*>::iterator poly = subTrias.begin() ; poly != subTrias.end(); ++poly)
  {
      domTemp = (*poly)->getDomainNumber();
      
      if( domTemp == 0 )
      {
        for(ii=0; ii<(*poly)->GetNumVertices(); ii++)
        {
          ptTemp = (*poly)->GetPoint(ii);
          
          dist = DistanceToLine2D(pt1, pt2, ptTemp);
          //cout << " dist = " << dist << endl;
          
          // if the distance is zero (within the tolerance) the the point 'ptTemp' is on the line
          if( abs(dist) <= tol )
          {
            //cout << " dist = " << dist << '\t' << " true " << endl;
            inserted = pointExists(myLine, ptTemp);

            if(!inserted)
            {
              myLine.push_back(ptTemp);
            }
          }
        }
      }
  }
*/

//

    //cout << " side = " << side << endl;

    myLine  lineTemp(pt1, pt2);
    vector<myPoint>  vecPts;

    lineTemp.computeNormal();
    lineTemp.computeAABB();

    domTemp = *std::max_element(domNums.begin(), domNums.end() ) - 1;

    vector<myPoint>  ptOut;

    //cout << " domTemp = " << domTemp << endl;
    //printf(" \t %12.8f \t %12.8f \t %12.8f \n ", pt1[0], pt1[1], pt1[2] );
    //printf(" \t %12.8f \t %12.8f \t %12.8f \n ", pt2[0], pt2[1], pt2[2] );
    //bbox.printSelf();
    //cout << GeomData->immSolidPtrs[domTemp]->ImmersedFaces.size() << endl;

    for(ii=0; ii<GeomData->immSolidPtrs[domTemp]->ImmersedFaces.size(); ii++)
    {
        // check if the boundingbox intersects any of the line segments
        // if it intersects then find the intersection points

        if( GeomData->immSolidPtrs[domTemp]->ImmersedFaces[ii]->doIntersect(bbox) )
        {
          //cout << " face " << ii << " intersects with cell AABB " << endl;
          GeomData->immSolidPtrs[domTemp]->ImmersedFaces[ii]->IntersectWithLine(lineTemp, ptOut);
        }
    } // for(ii=0; ii<ImmersedFaces.size(); ii++)


    int  d1 = GeomData->within(pt1);
    int  d2 = GeomData->within(pt2);

    //cout << " domTemp = " << domTemp << '\t' << d1 << '\t' << d2 << endl;

    BoundaryQuadrature[side].reset();

    if(d1 == d2)
      return  1;

    double  len;

    if(!d1 && d2) // pt1 in the fluid domain and pt2 in the domain 'domTemp'
    {
      vecPts.push_back(pt1);
      vecPts.push_back(ptOut[0]);

      ptTemp = ptOut[0] - pt1;
      len = ptTemp.norm();
    }
    if(d1 && !d2) // pt2 in the fluid domain and pt1 in the domain 'domTemp'
    {
      vecPts.push_back(ptOut[0]);
      vecPts.push_back(pt2);

      ptTemp = pt2 - ptOut[0];
      len = ptTemp.norm();
    }

    //cout << " side = " << side << endl;
    //for(ii=0; ii<2; ii++)
      //cout << vecPts[ii][0] << '\t' << vecPts[ii][1] << endl;

    //cout << " length = " << len << endl;

    switch(side)
    {
        case 0:
        case 1:

              nGauss = GeomData->gausspoints2.size();
              len = len/nGauss;

              BoundaryQuadrature[side].gausspoints.resize(nGauss);
              BoundaryQuadrature[side].gaussweights.resize(nGauss);

              for(ii=0; ii<nGauss; ii++)
              {
                // coordinates of Gauss points in physcial domain
                ptTemp = vecPts[0] + 0.5*(1.0 + GeomData->gausspoints2[ii])*(vecPts[1] - vecPts[0]);

                // physcial domain to parametric domain
                param = GeomData->computeParam(1, ptTemp[1]);

                // parametric domain to integration master-quadrilateral domain
                val = (2.0*param - knotSum[1])/knotIncr[1];

                BoundaryQuadrature[side].gausspoints[ii][0] = param;
                BoundaryQuadrature[side].gausspoints[ii][1] = val;

                BoundaryQuadrature[side].gaussweights[ii] = len*GeomData->gaussweights2[ii];
              }

        break;

        case 2:
        case 3:

              nGauss = GeomData->gausspoints1.size();
              len = len/nGauss;

              BoundaryQuadrature[side].gausspoints.resize(nGauss);
              BoundaryQuadrature[side].gaussweights.resize(nGauss);

              for(ii=0; ii<nGauss; ii++)
              {
                // coordinates of Gauss points in physcial domain
                ptTemp = vecPts[0] + 0.5*(1.0 + GeomData->gausspoints1[ii])*(vecPts[1] - vecPts[0]);

                // physcial domain to parametric domain
                param = GeomData->computeParam(0, ptTemp[0]);

                // parametric domain to integration master-quadrilateral domain
                val = (2.0*param - knotSum[0])/knotIncr[0];

                BoundaryQuadrature[side].gausspoints[ii][0] = val;
                BoundaryQuadrature[side].gausspoints[ii][1] = param;

                BoundaryQuadrature[side].gaussweights[ii] = len*GeomData->gaussweights1[ii];
              }

        break;

        default :

            cout << " Invalid 'side' value in GeomDataHBSplines::getBoundaryGPs2D " << endl;
        break;
    } //switch(side)
//
    //printVector(boundaryGPs1);
    //printVector(boundaryGPs2);

    //cout << " length = " << len << endl;

    return  1;
}


template<>
void TreeNode<1>::applyDirichletBCsCutFEMFluid(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
  return;
}


template<>
void TreeNode<2>::applyDirichletBCsCutFEMFluid(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
  if(DirichletData.size() > 0)
  {
    //if( isCutElement() )
    //{
      //vector<int>  domTemp;
      //std::vector<int>::iterator  itTemp;

      //for(vector<myPoly*>::iterator poly = subTrias.begin() ; poly != subTrias.end(); ++poly)
      //{
        //domTemp.push_back( (*poly)->getDomainNumber() );
      //}
      //printVector(domTemp);

      //if( ! any_of( domTemp.begin(), domTemp.end(), [](int aa){return aa==0;} ) )
      //{
        //cout << " true " << endl;
        //return;
      //}
    //}


    int ii, jj, aa, gp, nGauss, TI, TIp1, TIp2, index, TJ, TJp1, TJp2, dir, side;
    double  theta, y0, y1, Ta[3], Tb[3], res, JacTemp, temp, NitscheFact, pres;
    double  dvol, specVal, PENALTY, Jac, fact, rad, R, bb1, bb2, a;
    bool  isNitsche;

    bool   axsy = ((int)elmDat[2] == 1);
    double  rho = elmDat[3];
    double  mu  = elmDat[4];

    double  hx = bbox.maxBB[0]-bbox.minBB[0];
    double  hy = bbox.maxBB[1]-bbox.minBB[1];

    VectorXd  N(totnlbf), dN_dx(totnlbf), dN_dy(totnlbf);
    VectorXd  NN(totnlbf), dNN_dx(totnlbf), dNN_dy(totnlbf), trac(2);
    MatrixXd  grad(2,2), stress(2,2);
    myPoint  param, normal, geom;
    double *gws;
    myPoint *gps;
<<<<<<< HEAD:HBsplines/TreeNode9.cpp
    
    y0 = 0.0;    y1 = 0.41;
    y0 = -5.0;  y1 = 5.0;
    //y0 = 0.90;  y1 = 1.30;
    //y0 = 0.5;  y1 = 1.50;
    //y0 = 0.0;    y1 = 1.61;
=======
>>>>>>> collabchandan:src/HBsplines/TreeNode8.cpp

    //KimMoinFlow  analy(rho, mu);
    //Kovasznay  analy;
    //Stokes2DEx1  analy;

    for(aa=0;aa<DirichletData.size();aa++)
    {
        isNitsche   = false;
        side        = (int) (DirichletData[aa][0] - 1);
        dir         = (int) (DirichletData[aa][1] - 1);
        specVal     = DirichletData[aa][2];
        PENALTY     = DirichletData[aa][3];
        isNitsche   = ( (int) DirichletData[aa][4] == 1 );
        NitscheFact = DirichletData[aa][5];

        //PENALTY    = 1.0/max(hx, hy);                     // GeomData-> degree;

        //for symmetric Nitsche method   -> NitscheFact =  1.0
        //for unsymmetric Nitsche method -> NitscheFact = -1.0

        normal = GeomData->boundaryNormals[side];

        if(domNums.size() > 1)
        {
<<<<<<< HEAD:HBsplines/TreeNode9.cpp
          //TreeNode<2>::computeGaussPointsAdapIntegrationBoundary(side, 0, 0, 0, 0);

=======
>>>>>>> collabchandan:src/HBsplines/TreeNode8.cpp
          nGauss = BoundaryQuadrature[side].gausspoints.size() ;

          gps = &(BoundaryQuadrature[side].gausspoints[0]);
          gws = &(BoundaryQuadrature[side].gaussweights[0]);

          JacTemp = 1.0;
        }
        else
        {
          nGauss = GeomData->boundaryQuadrature2D[side].gausspoints.size();

          gps = &(GeomData->boundaryQuadrature2D[side].gausspoints[0]);
          gws = &(GeomData->boundaryQuadrature2D[side].gaussweights[0]);

          JacTemp = GeomData->boundaryJacobians[side][level];
        }

        for(gp=0; gp<nGauss; gp++)
        {
            param[0]  = 0.5*(knotIncr[0] * gps[gp][0] + knotSum[0]);
            param[1]  = 0.5*(knotIncr[1] * gps[gp][1] + knotSum[1]);

            dvol = JacTemp * gws[gp] ;

            geom[0] = GeomData->computeCoord(0, param[0]);
            geom[1] = GeomData->computeCoord(1, param[1]);

            if(axsy)
            {
              //if(side != 0)
                dvol *= 2.0*PI*geom[0];
            }

            if(parent == NULL)
            {
                GeomData->computeBasisFunctions2D(knotBegin, knotIncr, param, N, dN_dx, dN_dy );
            }
            else
            {
                GeomData->computeBasisFunctions2D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy );

                N = SubDivMat*NN;
                dN_dx = SubDivMat*dNN_dx;
                dN_dy = SubDivMat*dNN_dy;
            }

            specVal = DirichletData[aa][2];

              /*
              if(side == 0 )
              {
                if(dir == 0)
                {
                  // beam Bathe
                  //y0 = 0.0;    y1 = 0.5;
                  //specVal = DirichletData[aa][2]*(1.0-(50.0-geom[1])*(50.0-geom[1])/2500.0);

                  // vertical beam - Wall
                  //y0 = 0.0;    y1 = 0.6;
                  //specVal = DirichletData[aa][2]*(6.0/y1/y1)*(y1-geom[1])*(geom[1]-y0);

                  // Turek beam
                  y0 = 0.0;    y1 = 0.41;
                  specVal = DirichletData[aa][2]*(6.0/y1/y1)*(y1-geom[1])*(geom[1]-y0);

                  //specVal = DirichletData[aa][2]*(y1*y1-geom[1]*geom[1])/0.5625;
<<<<<<< HEAD:HBsplines/TreeNode9.cpp
                  //specVal = DirichletData[aa][2]*(6.0/(y1-y0)/(y1-y0))*(y1-geom[1])*(geom[1]-y0); // throttle valve
=======

                  // throttle valve
                  //specVal = DirichletData[aa][2]*(6.0/(y1-y0)/(y1-y0))*(y1-geom[1])*(geom[1]-y0);
>>>>>>> collabchandan:src/HBsplines/TreeNode8.cpp
                  //if(geom[1] >= y0 && geom[1] <= y1)

                  // heart-valve benchmark
                  //y0 = 0.0;    y1 = 1.61;
                  //specVal = DirichletData[aa][2]*(y1-geom[1])*(geom[1]-y0);

                  // single-leaf benchmark
                  //y0 = 0.0;    y1 = 2.0;
                  //specVal = DirichletData[aa][2]*(y1-geom[1])*(geom[1]-y0);

                  // heart-valve contacts
                  //y0 = 0.0;    y1 = 1.85;
                  //specVal = DirichletData[aa][2]*(y1-geom[1])*(geom[1]-y0);

                  // Neumann problem
                  //y0 = 0.0;    y1 = 1.0;
                  //specVal = DirichletData[aa][2]*(y1-geom[1])*(geom[1]-y0);

                  //else
                    //specVal = 0.0;
                    //specVal = DirichletData[aa][2]*(1.0-geom[1]*geom[1]);
<<<<<<< HEAD:HBsplines/TreeNode9.cpp
                }
              }
              */

              /*
              if(side == 3)
              {
                if(dir == 0)
                {
                  if(param[0]<0.5)
                    specVal = tanh(param[0]/knots[0][2]);
                  else
                    specVal = -tanh((param[0]-1.0)/knots[0][2]);
=======
>>>>>>> collabchandan:src/HBsplines/TreeNode8.cpp
                }
              }
              */

            //specVal = GeomData->analyDBC->computeValue(dir, xx, yy);
            //specVal = analy.computeValue(dir, geom[0], geom[1]);

            specVal *= timeFunction[0].prop;

            res = specVal - computeValue(dir, N);

            grad(0,0) = computeValue(0, dN_dx);
            grad(0,1) = computeValue(0, dN_dy);
            grad(1,0) = computeValue(1, dN_dx);
            grad(1,1) = computeValue(1, dN_dy);
            pres      = computeValue(2, N);

            stress = mu*grad;
            stress(0,0) -= pres;
            stress(1,1) -= pres;

            trac[0] = stress(0,0)*normal[0] + stress(0,1)*normal[1] ;
            trac[1] = stress(1,0)*normal[0] + stress(1,1)*normal[1] ;

            for(ii=0;ii<totnlbf2;ii++)
            {
                fact = (dvol*PENALTY)* N[ii] ;

                TI = ndof*ii+dir;

                Flocal(TI) += fact*res;

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

                    bb1 = N[ii]*dvol;

                    Ta[0] = (mu*dvol)*(normal[0]*dN_dx(ii)+normal[1]*dN_dy(ii));
                    Ta[1] = 0.0;
                    Ta[2] = -normal[0]*bb1;

                    for(jj=0;jj<totnlbf2;jj++)
                    {
                      TJ   = ndof*jj;
                      TJp1 = TJ+1;
                      TJp2 = TJ+2;

                      Tb[0] = mu*(normal[0]*dN_dx(jj)+normal[1]*dN_dy(jj));
                      Tb[1] = 0.0;
                      Tb[2] = -normal[0]*N(jj);

                      Klocal(TI, TJ)   -= (bb1*Tb[0]);
                      Klocal(TI, TJp1) -= (bb1*Tb[1]);
                      Klocal(TI, TJp2) -= (bb1*Tb[2]);

                      // Nitsche terms
                      Klocal(TI,   TJ) -= (Ta[0]*N(jj))*NitscheFact;
                      Klocal(TIp1, TJ) -= (Ta[1]*N(jj))*NitscheFact;
                      Klocal(TIp2, TJ) -= (Ta[2]*N(jj))*NitscheFact;
                    }

                    Flocal(TI)   += (bb1*trac[0]);

                    Flocal(TI)   -= (Ta[0]*res)*NitscheFact;
                    Flocal(TIp1) -= (Ta[1]*res)*NitscheFact;
                    Flocal(TIp2) -= (Ta[2]*res)*NitscheFact;
                  }
                }
                if(dir == 1)
                {
                  for(ii=0;ii<totnlbf2;ii++)
                  {
                    TI   = ndof*ii;
                    TIp1 = TI+1;
                    TIp2 = TI+2;

                    bb1 = N[ii]*dvol;

                    Ta[0] = 0.0;
                    Ta[1] = (mu*dvol)*(normal[0]*dN_dx(ii)+normal[1]*dN_dy(ii));
                    Ta[2] = -normal[1]*bb1;

                    for(jj=0;jj<totnlbf2;jj++)
                    {
                      TJ   = ndof*jj;
                      TJp1 = TJ+1;
                      TJp2 = TJ+2;

                      Tb[0] = 0.0;
                      Tb[1] = mu*(normal[0]*dN_dx(jj)+normal[1]*dN_dy(jj));
                      Tb[2] = -normal[1]*N(jj);

                      Klocal(TIp1, TJ)    -= (bb1*Tb[0]);
                      Klocal(TIp1, TJp1)  -= (bb1*Tb[1]);
                      Klocal(TIp1, TJp2)  -= (bb1*Tb[2]);

                      // Nitsche terms
                      Klocal(TI,   TJp1)  -= (Ta[0]*N(jj))*NitscheFact;
                      Klocal(TIp1, TJp1)  -= (Ta[1]*N(jj))*NitscheFact;
                      Klocal(TIp2, TJp1)  -= (Ta[2]*N(jj))*NitscheFact;
                    }

                    Flocal(TIp1) += (bb1*trac[1]);

                    Flocal(TI)   -= (Ta[0]*res)*NitscheFact;
                    Flocal(TIp1) -= (Ta[1]*res)*NitscheFact;
                    Flocal(TIp2) -= (Ta[2]*res)*NitscheFact;
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
void TreeNode<1>::applyNeumannBCsCutFEMFluid(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
  return;
}



template<>
void TreeNode<2>::applyNeumannBCsCutFEMFluid(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
    if( NeumannData.size() > 0 )
    {
        int ii, jj, aa, gp, nGauss, TI, TIp1, TIp2, index, dir, side;
        double  y0, y1, res, JacTemp;
        double  dvol, specVal, rad, freq, pres, PENALTY;

        VectorXd  NN(totnlbf), dNN_dx(totnlbf), dNN_dy(totnlbf);
        VectorXd  N(totnlbf2);

        myPoint  param, geom, normal;
        double *gws;
        myPoint *gps;

        bool   axsy = ((int)elmDat[2] == 1);
        double  rho = elmDat[3];
        double  mu  = elmDat[4];

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
                nGauss = GeomData->boundaryQuadrature2D[side].gausspoints.size();

                gps = &(GeomData->boundaryQuadrature2D[side].gausspoints[0]);
                gws = &(GeomData->boundaryQuadrature2D[side].gaussweights[0]);

                JacTemp = GeomData->boundaryJacobians[side][level];
            }

            for(gp=0; gp<nGauss; gp++)
            {
                param[0]  = 0.5*(knotIncr[0] * gps[gp][0] + knotSum[0]);
                param[1]  = 0.5*(knotIncr[1] * gps[gp][1] + knotSum[1]);

                dvol = JacTemp * gws[gp] ;

                geom[0] = GeomData->computeCoord(0, param[0]);
                geom[1] = GeomData->computeCoord(1, param[1]);

                if(axsy)
                {
                    if(side != 0)
                      dvol *= 2.0*PI*geom[0];
                }

                GeomData->computeBasisFunctions2D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy );

                if(parent == NULL)
                {
                  N = NN;
                }
                else
                {
                  N = SubDivMat*NN;
                }

                res = specVal * timeFunction[0].prop;

                //res = dvol * specVal * tanh(20.0*mpapTime.cur);
                //cout << specVal << '\t' << axsy << '\t' << dir << endl;
                //res *= (0.5*(1.0-cos(2.0*PI*SolnData->ElemProp.data[6]*mpapTime.cur)));
                //w1 = tanh(SolnData->ElemProp.data[5]*mpapTime.cur);

                //freq = 2.0*PI*10.0;
                //res = dvol* specVal * (0.5*(1.0-cos(freq*mpapTime.cur)));
                //res = dvol * specVal * sin(freq*mpapTime.cur);
                //res *= tanh(50.0*timeFunction[0].prop);

                res *= dvol;
                for(ii=0;ii<totnlbf2;ii++)
                  Flocal(ndof*ii+dir) += (res*N(ii));

            } // for(gp=0...
        } // for(aa=0;aa<NeumannData.size();aa++)
    } // if(NeumannData.size() > 0)

    return;
}



template<>
void TreeNode<1>::applyDerivativeBCsCutFEMFluid(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
  return;
}



template<>
void TreeNode<2>::applyDerivativeBCsCutFEMFluid(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
    if( DerivativeBCData.size() > 0 )
    {
        int ii, jj, aa, gp, nGauss, TI, TIp1, TIp2, index, dir, side;
        double  y0, y1, res, JacTemp;
        double  dvol, specVal, rad, freq, pres, PENALTY;

        VectorXd  N(totnlbf), dN_dx(totnlbf), dN_dy(totnlbf);
        VectorXd  NN(totnlbf), dNN_dx(totnlbf), dNN_dy(totnlbf), grad(2);
        MatrixXd  D(forAssyVec.size(),1);
        D.setZero();

        myPoint  param, geom, normal;
        double *gws;
        myPoint *gps;

        bool   axsy = ((int)elmDat[2] == 1);
        double  rho = elmDat[3];
        double  mu  = elmDat[4];

        for(aa=0;aa<DerivativeBCData.size();aa++)
        {
            side    = (int) (DerivativeBCData[aa][0] - 1);
            dir     = (int) (DerivativeBCData[aa][1] - 1);
            specVal = DerivativeBCData[aa][2];
            PENALTY = DerivativeBCData[aa][3];

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
                nGauss = GeomData->boundaryQuadrature2D[side].gausspoints.size();

                gps = &(GeomData->boundaryQuadrature2D[side].gausspoints[0]);
                gws = &(GeomData->boundaryQuadrature2D[side].gaussweights[0]);

                JacTemp = GeomData->boundaryJacobians[side][level];
            }


            for(gp=0; gp<nGauss; gp++)
            {
                param[0]  = 0.5*(knotIncr[0] * gps[gp][0] + knotSum[0]);
                param[1]  = 0.5*(knotIncr[1] * gps[gp][1] + knotSum[1]);

                dvol = JacTemp * gws[gp] ;

                geom[0] = GeomData->computeCoord(0, param[0]);
                geom[1] = GeomData->computeCoord(1, param[1]);

                if(axsy)
                {
                  if(side != 0)
                    dvol *= 2.0*PI*geom[0];
                }

                dvol *= PENALTY;

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

                if(dir == 0)
                {
                    grad(0) = computeValue(0, dN_dx);
                    grad(1) = computeValue(0, dN_dy);

                    for(ii=0;ii<totnlbf2;ii++)
                    {
                        TI = ndof*ii;
                        D(TI)   = (dN_dx[ii]*normal[0] + dN_dy[ii]*normal[1]);
                        D(TI+1) = 0.0;
                    }

                    res = specVal - (grad(0)*normal[0]+grad(1)*normal[1]);
                }

                if(dir == 1)
                {
                    grad(0) = computeValue(1, dN_dx);
                    grad(1) = computeValue(1, dN_dy);

                    for(ii=0;ii<totnlbf2;ii++)
                    {
                        TI = ndof*ii;
                        D(TI)   = 0.0;
                        D(TI+1) = (dN_dx[ii]*normal[0] + dN_dy[ii]*normal[1]);
                    }

                    res = specVal - (grad(0)*normal[0]+grad(1)*normal[1]);
                }

                Flocal += (dvol*res)*D;
                Klocal += (dvol*D)*D.transpose();

            } // for(gp=0...
        } // for(aa=0;aa<NeumannData.size();aa++)
    } // if(NeumannData.size() > 0)

    return;
}



