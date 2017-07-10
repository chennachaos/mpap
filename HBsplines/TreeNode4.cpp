
#include "TreeNode.h"
#include "MpapTime.h"
#include "Functions.h"
#include "DistFunctions.h"
#include "TimeFunction.h"
#include "SolutionData.h"
#include "myDataIntegrateCutFEM.h"

#include "BasisFunctionsBSpline.h"


extern MpapTime mpapTime;
extern List<TimeFunction> timeFunction;



/*
template<>
void TreeNode<2>::calcStiffnessAndResidualGFEM(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
    // GFEM for Navier-Stokes flow 
    // full subroutine
    //
    //////////////////////////////////////////////////
    
    Circle  circle2(0.0, 0.0, 0.45);

    int ii, jj, gp1, gp2, TI, TIp1, TIp2, count, TJ, TJp1, TJp2, nU, nP;
   
    double  uu, vv, Jac, dvol, fact, b1, b2, b3, b4, b5, b6, b7, b8, xx, yy, acceFact, HH, dist, PERM;
    double  pres, Da, Db, af, am, d1, c1, muTaf, rad, urdr, urdr2, nu, alpha1, alpha2, conv, tau[2];

    nU  = forAssyVec.size();
    nP  = forAssyVec2.size();
    count = nU + nP;

    VectorXd  NN(totnlbf), dNN_dx(totnlbf), d2NN_dx2(totnlbf), dNN_dy(totnlbf), d2NN_dy2(totnlbf), vectmp(totnlbf);
    VectorXd  N, dN_dx, d2N_dx2, dN_dy, d2N_dy2;
    VectorXd  res(3), dp(2), Du(2), vel(2), velDot(2), force(2);
    MatrixXd  F(2,2), FN(2,2), stress(2,2);

    af = SolnData->td(2);
    am = SolnData->td(3);
    acceFact = rho*am*SolnData->td(9);
    muTaf = mu*af;
    
    nu = 0.4999;
    nu = 0.5;
    b1 = 1.5*(1.0-2.0*nu)/(1.0+nu);
    alpha1 = 1.0 - b1/3.0;
    //alpha1 = 1.0;
    alpha2 = b1/mu;
    
    conv = 1.0;
    
    tau[0] = stabParam * 1.0;

    count = 0;
    for(gp2=0;gp2<GeomData->getNGP(1);gp2++)
    {
       vv  = 0.5*(knots[1][2] * GeomData->gausspoints2[gp2] + knots[1][3]);
       Jac = GeomData->gaussweights2[gp2] * JacMultElem;
       
       for(gp1=0;gp1<GeomData->getNGP(0);gp1++)
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

          xx = GeomData->computeCoord(0, uu);
          yy = GeomData->computeCoord(1, vv);
          
          dist = circle2.ComputeDistance(xx, yy);
          if(dist > 0.0)
            PERM = 0.0;
          else
            PERM = 1.0e2;
          
          PERM = 0.0;

          //muTaf = (mu + PERM)*af;

          //cout << " dist = " << dist << endl;

          //if(circle2.checkPointLocation(xx, yy))
            //stabParam2 = 0.0;
          //else
            //stabParam2 = stabParam ;

          vel(0) = computeValueCur(0, N);
          vel(1) = computeValueCur(1, N);

          F(0,0) = computeValueCur(0, dN_dx);
          F(0,1) = computeValueCur(0, dN_dy);
          F(1,0) = computeValueCur(1, dN_dx);
          F(1,1) = computeValueCur(1, dN_dy);
          
          pres   = computeValue(2, N);
          dp(0)  = computeValue(2, dN_dx);
          dp(1)  = computeValue(2, dN_dy);

          //pres   = computeValue2(0, N);
          //dp(0)  = computeValue2(0, dN_dx);
          //dp(1)  = computeValue2(0, dN_dy);

          //vectmp = d2N_dx2 + d2N_dy2;
          //Du(0) = computeValue(0, vectmp);
          //Du(1) = computeValue(1, vectmp);

          velDot(0) = computeValueDotCur(0, N);
          velDot(1) = computeValueDotCur(1, N);
          
          // this is pseudo-stress
          //stress = (mu + PERM)*F;
          stress = mu *F;
          stress(0,0) -= pres*alpha1;
          stress(1,1) -= pres*alpha1;

          force.setZero();
          //force(0) = analy.computeForce(0, xx, yy);
          //force(1) = analy.computeForce(1, xx, yy);
          //cout << force(0) << '\t' << force(1) << endl;

          //force(0) = computeForce(0, N);
          //force(1) = computeForce(1, N);

          velDot += (F*vel)*conv;
          velDot *= rho;

          if(axsy)
          {
            rad = xx;

            urdr  = vel(0)/rad;
            urdr2 = urdr/rad;
            dvol *= (2.0*PI*rad);
          }

        dist = 1.0;
        if(dist > 0.0) // fluid
        {
          //muTaf = (mu + PERM)*af;
          for(ii=0;ii<totnlbf2;ii++)
          {
             TI   = ndof*ii;
             TIp1 = TI+1;
             //TIp2 = nU+ii;
             TIp2 = TI+2;

             b1 = dN_dx[ii]*dvol;
             b2 = dN_dy[ii]*dvol;
             b4 = N[ii]*dvol;

             b5 = muTaf*b1;
             b6 = muTaf*b2;
             b8 = af*b4;

             c1 = acceFact*b4;

             FN = (rho*b8)*F*conv;

             for(jj=0;jj<totnlbf2;jj++)
             {
               TJ   = ndof*jj;
               TJp1 = TJ+1;
               //TJp2 = nU+jj;
               TJp2 = TJ+2;

               // time acceleration term
               fact = c1*N(jj) ;

               // diffusion term
               fact += b5*dN_dx(jj)+b6*dN_dy(jj);
               
               // Brinkmann flow term
               fact += PERM * b8*N(jj);

               Klocal(TI,   TJ)   += fact;
               Klocal(TIp1, TJp1) += fact;

               // convection term

               Db = rho*(vel(0)*dN_dx(jj) + vel(1)*dN_dy(jj))*conv;

               Klocal(TI,   TJ)   += (FN(0,0)*N(jj) + b8*Db);
               Klocal(TI,   TJp1) += (FN(0,1)*N(jj));
               Klocal(TIp1, TJ)   += (FN(1,0)*N(jj));
               Klocal(TIp1, TJp1) += (FN(1,1)*N(jj) + b8*Db);

               // pressure term
               Klocal(TI,   TJp2) -= (b1*N(jj)*alpha1);
               Klocal(TIp1, TJp2) -= (b2*N(jj)*alpha1);

               // continuity equation
               Klocal(TIp2, TJ)   -= (b8*dN_dx(jj));
               Klocal(TIp2, TJp1) -= (b8*dN_dy(jj));
               Klocal(TIp2, TJp2) += (b8*N(jj)*alpha2);

               //stabilisation terms
               Klocal(TIp2, TJp2) -= (b1*dN_dx(jj)+b2*dN_dy(jj))*tau[0];

               if(axsy)
               {
                  // diffusion term
                  Klocal(TI, TJ)     += (mu*b4*N(jj)/rad/rad);
                  Klocal(TI, TJp2)   -= (b4*N(jj)/rad);

                  // continuity equation
                  Klocal(TIp2, TJ)   -= (b4*N(jj)/rad);
               }
             }

             Flocal(TI)   += (b4*(force(0)-velDot(0)-PERM*vel(0)) - b1*stress(0,0) - b2*stress(0,1) );
             Flocal(TIp1) += (b4*(force(1)-velDot(1)-PERM*vel(1)) - b1*stress(1,0) - b2*stress(1,1) );
             Flocal(TIp2) += (b4*(F.trace()-pres*alpha2));

             // stabilisation terms
             Flocal(TIp2) += (tau[0]*(b1*dp(0)+b2*dp(1)));

             if(axsy)
             {
                Flocal(TI)   += (-b4*(mu*vel(0)/rad/rad));
                Flocal(TI)   += (b4*pres/rad);
                Flocal(TIp2) += (b4*vel(0)/rad);
             }
          }
        }
        else
        {
          for(ii=0;ii<totnlbf2;ii++)
          {
             TI   = ndof*ii;
             TIp1 = TI+1;
             //TIp2 = nU+ii;
             TIp2 = TI+2;

             b4 = N[ii]*dvol;

             b8 = af*b4;

             for(jj=0;jj<totnlbf2;jj++)
             {
               TJ   = ndof*jj;
               TJp1 = TJ+1;
               //TJp2 = nU+jj;
               TJp2 = TJ+2;

               fact = b8*N(jj);

               Klocal(TI,   TJ)   += fact;
               Klocal(TIp1, TJp1) += fact;
               //Klocal(TIp2, TJp2) += fact;

               // continuity equation
               //Klocal(TIp2, TJ)   += (b8*dN_dx(jj));
               //Klocal(TIp2, TJp1) += (b8*dN_dy(jj));

             }

             Flocal(TI)   += (b4*(0.0-vel(0)));
             Flocal(TIp1) += (b4*(0.0-vel(1)));
             //Flocal(TIp2) += (b4*(0.0-pres));
             //Flocal(TIp2) -= (b4*(F.trace()));
          }
        }
    }//gp1
    }//gp2
    
    return;
}
*/


/*
template<>
void TreeNode<2>::calcStiffnessAndResidualGFEM(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
    // GFEM for Navier-Stokes flow 
    // without any stabilisation
    // fullly implicit
    //////////////////////////////////////////////////

    int ii, jj, gp1, gp2, TI, TIp1, TIp2, count, TJ, TJp1, TJp2;
   
    double  Jac, dvol, fact, b1, b2, b3, b4, b5, b6, b7, b8;
    double  pres, Da, Db, af, am, d1, c1, muTaf, rad, urdr, urdr2, acceFact;

    VectorXd  NN(totnlbf), dNN_dx(totnlbf), dNN_dy(totnlbf);
    VectorXd  N, dN_dx, dN_dy;
    VectorXd  res(3), dp(2), Du(2), vel(2), velDot(2), force(2), res2(2), gradTvel(2);
    MatrixXd  grad(2,2), gradN(2,2), stress(2,2);
    myPoint  param, geom;

    bool   axsy = ((int)elmDat[2] == 1);
    double  rho = elmDat[3];
    double  mu  = elmDat[4];

    af = SolnData->td(2);
    am = SolnData->td(3);
    acceFact = am*SolnData->td(9);
    muTaf = mu*af;

    for(gp2=0;gp2<GeomData->getNGP(1);gp2++)
    {
       param[1]  = 0.5*(knots[1][2] * GeomData->gausspoints2[gp2] + knots[1][3]);
       Jac = GeomData->gaussweights2[gp2] * JacMultElem;
       
       for(gp1=0;gp1<GeomData->getNGP(0);gp1++)
       {
          param[0]   = 0.5*(knots[0][2] * GeomData->gausspoints1[gp1] + knots[0][3]);
          dvol = GeomData->gaussweights1[gp1] * Jac;

          //GeomData->computeBasisFunctions2D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy, d2NN_dx2, d2NN_dy2);
          
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

          geom[0] = GeomData->computeCoord(0, param[0]);
          geom[1] = GeomData->computeCoord(1, param[1]);

          vel(0) = computeValueCur(0, N);
          vel(1) = computeValueCur(1, N);

          grad(0,0) = computeValueCur(0, dN_dx);
          grad(0,1) = computeValueCur(0, dN_dy);
          grad(1,0) = computeValueCur(1, dN_dx);
          grad(1,1) = computeValueCur(1, dN_dy);

          pres   = computeValue(2, N);
          dp(0)  = computeValue(2, dN_dx);
          dp(1)  = computeValue(2, dN_dy);

          //pres   = computeValue2(0, N);
          //dp(0)  = computeValue2(0, dN_dx);
          //dp(1)  = computeValue2(0, dN_dy);

          //vectmp = d2N_dx2 + d2N_dy2;
          //Du(0) = computeValue(0, vectmp);
          //Du(1) = computeValue(1, vectmp);

          velDot(0) = computeValueDotCur(0, N);
          velDot(1) = computeValueDotCur(1, N);
          
          // this is pseudo-stress
          stress = mu*grad;
          stress(0,0) -= pres;
          stress(1,1) -= pres;

          force.setZero();
          //force(0) = analy.computeForce(0, xx, yy);
          //force(1) = analy.computeForce(1, xx, yy);
          //cout << force(0) << '\t' << force(1) << endl;

          //force(0) = computeForce(0, N);
          //force(1) = computeForce(1, N);

          gradTvel = grad*vel ;

          res2(0) = rho*(velDot(0) + gradTvel(0) - force(0)) ;
          res2(1) = rho*(velDot(1) + gradTvel(1) - force(1)) ;

          if(axsy)
          {
            rad = geom[0];

            urdr  = vel(0)/rad;
            urdr2 = urdr/rad;
            dvol *= (2.0*PI*rad);
          }

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

             c1 = (rho*acceFact)*b4;

             for(jj=0;jj<totnlbf2;jj++)
             {
               TJ   = ndof*jj;
               TJp1 = TJ+1;
               TJp2 = TJ+2;

               // time acceleration term
               fact = c1*N(jj) ;

               // diffusion term
               fact += ( b5*dN_dx(jj)+b6*dN_dy(jj) );
               
               Klocal(TI,   TJ)   += fact;
               Klocal(TIp1, TJp1) += fact;

               // convection term

               Db = rho*(vel(0)*dN_dx(jj) + vel(1)*dN_dy(jj));

               gradN = grad*(rho*N(jj));

               gradN(0,0) += Db;
               gradN(1,1) += Db;

               Klocal(TI,   TJ)   += ( b8*gradN(0,0) );
               Klocal(TI,   TJp1) += ( b8*gradN(0,1) );
               Klocal(TIp1, TJ)   += ( b8*gradN(1,0) );
               Klocal(TIp1, TJp1) += ( b8*gradN(1,1) );

               // pressure term
               Klocal(TI,   TJp2) -= (b1*N(jj));
               Klocal(TIp1, TJp2) -= (b2*N(jj));

               // continuity equation
               Klocal(TIp2, TJ)   -= (b8*dN_dx(jj));
               Klocal(TIp2, TJp1) -= (b8*dN_dy(jj));

               if(axsy)
               {
                  // diffusion term
                  Klocal(TI, TJ)     += (b4 * (mu/rad/rad) * (af*N(jj)) );
                  Klocal(TI, TJp2)   -= (b4 * N(jj)/rad);

                  // continuity equation
                  Klocal(TIp2, TJ)   -= (b4 * af*N(jj)/rad);
               }
             }

             Flocal(TI)   -= (b4*res2(0) + b1*stress(0,0) + b2*stress(0,1) );
             Flocal(TIp1) -= (b4*res2(1) + b1*stress(1,0) + b2*stress(1,1) );
             Flocal(TIp2) += (b4*grad.trace());

             if(axsy)
             {
                Flocal(TI)   -= (b4 * (mu/rad/rad) * vel(0) );
                Flocal(TI)   += (b4 * pres/rad);
                Flocal(TIp2) += (b4 * vel(0)/rad);
             }
          }
    }//gp1
    }//gp2
    
    return;
}
*/



/*
template<>
void TreeNode<2>::calcStiffnessAndResidualGFEM(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
    // GFEM for Navier-Stokes flow 
    // without any stabilisation
    // semi-implicit  type A
    //////////////////////////////////////////////////

    int ii, jj, gp1, gp2, TI, TIp1, TIp2, count, TJ, TJp1, TJp2;
   
    double  Jac, dvol, fact, b1, b2, b3, b4, b5, b6, b7, b8;
    double  pres, Da, Db, af, am, d1, c1, muTaf, rad, urdr, urdr2, acceFact;

    VectorXd  NN(totnlbf), dNN_dx(totnlbf), dNN_dy(totnlbf);
    VectorXd  N, dN_dx, dN_dy;
    VectorXd  vel(2), velDot(2), force(2), res2(2), gradTvel(2), velPrev(2);
    MatrixXd  grad(2,2), gradPrev(2,2), gradN(2,2), stress(2,2);
    myPoint  param, geom;

    bool   axsy = ((int)elmDat[2] == 1);
    double  rho = elmDat[3];
    double  mu  = elmDat[4];

    af = SolnData->td(2);
    am = SolnData->td(3);
    acceFact = am*SolnData->td(9);
    muTaf = mu*af;

    for(gp2=0;gp2<GeomData->getNGP(1);gp2++)
    {
       param[1]  = 0.5*(knots[1][2] * GeomData->gausspoints2[gp2] + knots[1][3]);
       Jac = GeomData->gaussweights2[gp2] * JacMultElem;
       
       for(gp1=0;gp1<GeomData->getNGP(0);gp1++)
       {
          param[0]   = 0.5*(knots[0][2] * GeomData->gausspoints1[gp1] + knots[0][3]);
          dvol = GeomData->gaussweights1[gp1] * Jac;

          //GeomData->computeBasisFunctions2D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy, d2NN_dx2, d2NN_dy2);
          
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

          //geom[0] = GeomData->computeCoord(0, param[0]);
          //geom[1] = GeomData->computeCoord(1, param[1]);

          vel(0) = computeValueCur(0, N);
          vel(1) = computeValueCur(1, N);

          grad(0,0) = computeValueCur(0, dN_dx);
          grad(0,1) = computeValueCur(0, dN_dy);
          grad(1,0) = computeValueCur(1, dN_dx);
          grad(1,1) = computeValueCur(1, dN_dy);

          velPrev(0) = computeValuePrev(0, N);
          velPrev(1) = computeValuePrev(1, N);

          gradPrev(0,0) = computeValuePrev(0, dN_dx);
          gradPrev(0,1) = computeValuePrev(0, dN_dy);
          gradPrev(1,0) = computeValuePrev(1, dN_dx);
          gradPrev(1,1) = computeValuePrev(1, dN_dy);

          pres   = computeValue(2, N);

          velDot(0) = computeValueDotCur(0, N);
          velDot(1) = computeValueDotCur(1, N);
          
          // this is pseudo-stress
          stress = mu*grad;
          stress(0,0) -= pres;
          stress(1,1) -= pres;

          force.setZero();
          //force(0) = analy.computeForce(0, xx, yy);
          //force(1) = analy.computeForce(1, xx, yy);
          //cout << force(0) << '\t' << force(1) << endl;

          //force(0) = computeForce(0, N);
          //force(1) = computeForce(1, N);

          gradTvel = grad*velPrev;

          res2(0) = rho*(velDot(0) + gradTvel(0) - force(0)) ;
          res2(1) = rho*(velDot(1) + gradTvel(1) - force(1)) ;

          if(axsy)
          {
            rad = geom[0];

            urdr  = vel(0)/rad;
            urdr2 = urdr/rad;
            dvol *= (2.0*PI*rad);
          }

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

             c1 = (rho*acceFact)*b4;

             for(jj=0;jj<totnlbf2;jj++)
             {
               TJ   = ndof*jj;
               TJp1 = TJ+1;
               TJp2 = TJ+2;

               // time acceleration term
               fact = c1*N(jj) ;

               // diffusion term
               fact += ( b5*dN_dx(jj)+b6*dN_dy(jj) );
               
               Klocal(TI,   TJ)   += fact;
               Klocal(TIp1, TJp1) += fact;

              // convection term - semi-implicit type A

              Db = rho*(velPrev(0)*dN_dx(jj) + velPrev(1)*dN_dy(jj));

              Klocal(TI,   TJ)   += (b8*Db);
              Klocal(TIp1, TJp1) += (b8*Db);

               // pressure term
               Klocal(TI,   TJp2) -= (b1*N(jj));
               Klocal(TIp1, TJp2) -= (b2*N(jj));

               // continuity equation
               Klocal(TIp2, TJ)   -= (b8*dN_dx(jj));
               Klocal(TIp2, TJp1) -= (b8*dN_dy(jj));

               if(axsy)
               {
                  // diffusion term
                  Klocal(TI, TJ)     += (b4 * (mu/rad/rad) * (af*N(jj)) );
                  Klocal(TI, TJp2)   -= (b4 * N(jj)/rad);

                  // continuity equation
                  Klocal(TIp2, TJ)   -= (b4 * af*N(jj)/rad);
               }
             }

             Flocal(TI)   -= (b4*res2(0) + b1*stress(0,0) + b2*stress(0,1) );
             Flocal(TIp1) -= (b4*res2(1) + b1*stress(1,0) + b2*stress(1,1) );
             Flocal(TIp2) += (b4*grad.trace());

             if(axsy)
             {
                Flocal(TI)   -= (b4 * (mu/rad/rad) * vel(0) );
                Flocal(TI)   += (b4 * pres/rad);
                Flocal(TIp2) += (b4 * vel(0)/rad);
             }
          }
    }//gp1
    }//gp2
    
    return;
}
*/



//
template<>
void TreeNode<2>::calcStiffnessAndResidualGFEM(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
    // GFEM for Navier-Stokes flow 
    // without any stabilisation
    // semi-implicit  type B
    //////////////////////////////////////////////////

    int ii, jj, gp1, gp2, TI, TIp1, TIp2, count, TJ, TJp1, TJp2;
   
    double  Jac, dvol, fact, b1, b2, b3, b4, b5, b6, b7, b8;
    double  pres, Da, Db, af, am, d1, c1, muTaf, rad, urdr, urdr2, acceFact;

    VectorXd  NN(totnlbf), dNN_dx(totnlbf), dNN_dy(totnlbf);
    VectorXd  N, dN_dx, dN_dy;
    VectorXd  vel(2), velDot(2), force(2), res2(2), gradTvel(2), velPrev(2);
    MatrixXd  grad(2,2), gradPrev(2,2), gradN(2,2), stress(2,2);
    myPoint  param, geom;

    bool   axsy = ((int)elmDat[2] == 1);
    double  rho = elmDat[3];
    double  mu  = elmDat[4];

    af = SolnData->td(2);
    am = SolnData->td(3);
    acceFact = am*SolnData->td(9);
    muTaf = mu*af;

    for(gp2=0;gp2<GeomData->getNGP(1);gp2++)
    {
       param[1]  = 0.5*(knots[1][2] * GeomData->gausspoints2[gp2] + knots[1][3]);
       Jac = GeomData->gaussweights2[gp2] * JacMultElem;
       
       for(gp1=0;gp1<GeomData->getNGP(0);gp1++)
       {
          param[0]   = 0.5*(knots[0][2] * GeomData->gausspoints1[gp1] + knots[0][3]);
          dvol = GeomData->gaussweights1[gp1] * Jac;

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

          //geom[0] = GeomData->computeCoord(0, param[0]);
          //geom[1] = GeomData->computeCoord(1, param[1]);

          vel(0) = computeValueCur(0, N);
          vel(1) = computeValueCur(1, N);

          grad(0,0) = computeValueCur(0, dN_dx);
          grad(0,1) = computeValueCur(0, dN_dy);
          grad(1,0) = computeValueCur(1, dN_dx);
          grad(1,1) = computeValueCur(1, dN_dy);

          velPrev(0) = computeValuePrev(0, N);
          velPrev(1) = computeValuePrev(1, N);

          gradPrev(0,0) = computeValuePrev(0, dN_dx);
          gradPrev(0,1) = computeValuePrev(0, dN_dy);
          gradPrev(1,0) = computeValuePrev(1, dN_dx);
          gradPrev(1,1) = computeValuePrev(1, dN_dy);

          pres   = computeValue(2, N);

          velDot(0) = computeValueDotCur(0, N);
          velDot(1) = computeValueDotCur(1, N);
          
          // this is pseudo-stress
          stress = mu*grad;
          stress(0,0) -= pres;
          stress(1,1) -= pres;

          force.setZero();
          //force(0) = analy.computeForce(0, xx, yy);
          //force(1) = analy.computeForce(1, xx, yy);
          //cout << force(0) << '\t' << force(1) << endl;

          //force(0) = computeForce(0, N);
          //force(1) = computeForce(1, N);

          gradTvel = gradPrev*vel + grad*velPrev - gradPrev*velPrev;

          res2(0) = rho*(velDot(0) + gradTvel(0) - force(0)) ;
          res2(1) = rho*(velDot(1) + gradTvel(1) - force(1)) ;

          if(axsy)
          {
            rad = geom[0];

            urdr  = vel(0)/rad;
            urdr2 = urdr/rad;
            dvol *= (2.0*PI*rad);
          }

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

             c1 = (rho*acceFact)*b4;

             for(jj=0;jj<totnlbf2;jj++)
             {
               TJ   = ndof*jj;
               TJp1 = TJ+1;
               TJp2 = TJ+2;

               // time acceleration term
               fact = c1*N(jj) ;

               // diffusion term
               fact += ( b5*dN_dx(jj)+b6*dN_dy(jj) );
               
               Klocal(TI,   TJ)   += fact;
               Klocal(TIp1, TJp1) += fact;

               // convection term

               Db = rho*(velPrev(0)*dN_dx(jj) + velPrev(1)*dN_dy(jj));

               gradN = gradPrev*(rho*N(jj));

               gradN(0,0) += Db;
               gradN(1,1) += Db;

               Klocal(TI,   TJ)   += ( b8*gradN(0,0) );
               Klocal(TI,   TJp1) += ( b8*gradN(0,1) );
               Klocal(TIp1, TJ)   += ( b8*gradN(1,0) );
               Klocal(TIp1, TJp1) += ( b8*gradN(1,1) );

               // pressure term
               Klocal(TI,   TJp2) -= (b1*N(jj));
               Klocal(TIp1, TJp2) -= (b2*N(jj));

               // continuity equation
               Klocal(TIp2, TJ)   -= (b8*dN_dx(jj));
               Klocal(TIp2, TJp1) -= (b8*dN_dy(jj));

               if(axsy)
               {
                  // diffusion term
                  Klocal(TI, TJ)     += (b4 * (mu/rad/rad) * (af*N(jj)) );
                  Klocal(TI, TJp2)   -= (b4 * N(jj)/rad);

                  // continuity equation
                  Klocal(TIp2, TJ)   -= (b4 * af*N(jj)/rad);
               }
             }

             Flocal(TI)   -= (b4*res2(0) + b1*stress(0,0) + b2*stress(0,1) );
             Flocal(TIp1) -= (b4*res2(1) + b1*stress(1,0) + b2*stress(1,1) );
             Flocal(TIp2) += (b4*grad.trace());

             if(axsy)
             {
                Flocal(TI)   -= (b4 * (mu/rad/rad) * vel(0) );
                Flocal(TI)   += (b4 * pres/rad);
                Flocal(TIp2) += (b4 * vel(0)/rad);
             }
          }
    }//gp1
    }//gp2
    
    return;
}
//



/*
template<>
void TreeNode<2>::calcStiffnessAndResidualGFEM(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
    // GFEM for Navier-Stokes flow 
    // full subroutine
    // with proper stabilisation
    //////////////////////////////////////////////////


    int ii, jj, gp1, gp2, TI, TIp1, TIp2, count, TJ, TJp1, TJp2;
   
    double  Jac, dvol, fact, fact2, b1, b2, b3, b4, b5, b6, b7, b8, xx, yy, acceFact;
    double  pres, Da, Db, af, am, muTaf, rad, urdr, urdr2, h2, tau[3], eps;

    VectorXd  NN(totnlbf), dNN_dx(totnlbf), d2NN_dx2(totnlbf), dNN_dy(totnlbf), d2NN_dy2(totnlbf), d2N(totnlbf);
    VectorXd  N, dN_dx, d2N_dx2, dN_dy, d2N_dy2;
    VectorXd  res(3), res2(2), dp(2), Du(2), vel(2), velDot(2), force(2), Fvel(2), rStab(3);
    MatrixXd  Dj(2, 3), F(2,2), FN(2,2), stress(2,2);
    myPoint  param;
    Dj.setZero();

    bool   axsy = ((int)elmDat[2] == 1);
    double  rho = elmDat[3];
    double  mu  = elmDat[4];

    volume = GeomData->getGridLength(0) * knots[0][2] * GeomData->getGridLength(1) * knots[1][2];
//    cout << volume << '\t' << volume << endl;
    h2 = 4.0*volume/PI;

    double  stabParam = h2/(12.0*mu)/degree[0]/degree[1];
    tau[0] = elmDat[8]*stabParam;      // SUPG
    tau[1] = elmDat[9]*stabParam;//rho;  // PSPG
    tau[2] = elmDat[10]*stabParam*rho; // LSIC

    af = SolnData->td(2);
    am = SolnData->td(1);
    acceFact = am*SolnData->td(9);
    muTaf = mu*af;
    
    //Stokes2DEx2 analy;
    //Kovasznay  analy;
    //analy.SetPressure(0.0);

    count = 0;
    for(gp2=0;gp2<GeomData->getNGP(1);gp2++)
    {
      param[1] = 0.5*(knots[1][2] * GeomData->gausspoints2[gp2] + knots[1][3]);
      Jac = GeomData->gaussweights2[gp2] * JacMultElem;

      for(gp1=0;gp1<GeomData->getNGP(0);gp1++)
      {
          param[0]  = 0.5*(knots[0][2] * GeomData->gausspoints1[gp1] + knots[0][3]);
          dvol = GeomData->gaussweights1[gp1] * Jac;

          GeomData->computeBasisFunctions2D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy, d2NN_dx2, d2NN_dy2);

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

          xx = GeomData->computeCoord(0, param[0]);
          yy = GeomData->computeCoord(1, param[1]);

          //if(circle2.checkPointLocation(xx, yy))
            //stabParam2 = 0.0;
          //else
            //stabParam2 = stabParam ;

          vel(0) = computeValueCur(0, N);
          vel(1) = computeValueCur(1, N);

          F(0,0) = computeValueCur(0, dN_dx);
          F(0,1) = computeValueCur(0, dN_dy);
          F(1,0) = computeValueCur(1, dN_dx);
          F(1,1) = computeValueCur(1, dN_dy);
          
          pres   = computeValue(2, N);
          dp(0)  = computeValue(2, dN_dx);
          dp(1)  = computeValue(2, dN_dy);

          //pres   = computeValue2(0, N);
          //dp(0)  = computeValue2(0, dN_dx);
          //dp(1)  = computeValue2(0, dN_dy);

          d2N = d2N_dx2 + d2N_dy2;
          Du(0) = computeValueCur(0, d2N);
          Du(1) = computeValueCur(1, d2N);

          velDot(0) = computeValueDotCur(0, N);
          velDot(1) = computeValueDotCur(1, N);
          
          // this is pseudo-stress
          stress = mu*F;
          stress(0,0) -= pres;
          stress(1,1) -= pres;

          force.setZero();
          //force(0) = analy.computeForce(0, xx, yy);
          //force(1) = analy.computeForce(1, xx, yy);
          //cout << force(0) << '\t' << force(1) << endl;

          //force(0) = computeForce(0, N);
          //force(1) = computeForce(1, N);

          Fvel = F*vel;

          if(axsy)
          {
            rad = xx;

            urdr  = vel(0)/rad;
            urdr2 = urdr/rad;
            dvol *= (2.0*PI*rad);
          }

          res2(0) = rho*(velDot(0) + Fvel(0) - force(0)) ;
          res2(1) = rho*(velDot(1) + Fvel(1) - force(1)) ;

          rStab(0) = res2(0) - mu*Du(0) + dp(0) ;
          rStab(1) = res2(1) - mu*Du(1) + dp(1) ;

          //if( abs(xx) <= 0.0005)
            //stabParam2 = stabParam;
          //else
            //stabParam2 = 0.0;

          for(ii=0;ii<totnlbf2;ii++)
          {
             TI   = ndof*ii;
             TIp1 = TI+1;
             //TIp2 = nU+ii;
             TIp2 = TI+2;

             b1 = dN_dx[ii]*dvol;
             b2 = dN_dy[ii]*dvol;
             b4 = N[ii]*dvol;

             b5 = muTaf*b1;
             b6 = muTaf*b2;
             b8 = af*b4;
             
             Da = (vel(0)*b1 + vel(1)*b2)*tau[0];

             for(jj=0;jj<totnlbf2;jj++)
             {
               TJ   = ndof*jj;
               TJp1 = TJ+1;
               //TJp2 = nU+jj;
               TJp2 = TJ+2;

               fact2 = rho*acceFact*N(jj);
               
               // time acceleration term
               fact = b4*fact2 ;

               // diffusion term
               fact += b5*dN_dx(jj)+b6*dN_dy(jj);

               Klocal(TI,   TJ)   += fact;
               Klocal(TIp1, TJp1) += fact;

               // convection term

               FN = F*(rho*N(jj));

               Db = rho*(vel(0)*dN_dx(jj) + vel(1)*dN_dy(jj));
               
               FN(0,0) += Db;
               FN(1,1) += Db;

               Klocal(TI,   TJ)   += (b8*FN(0,0));
               Klocal(TI,   TJp1) += (b8*FN(0,1));
               Klocal(TIp1, TJ)   += (b8*FN(1,0));
               Klocal(TIp1, TJp1) += (b8*FN(1,1));

               // pressure term
               Klocal(TI,   TJp2) -= (b1*N(jj));
               Klocal(TIp1, TJp2) -= (b2*N(jj));

               // continuity equation
               Klocal(TIp2, TJ)   += (b8*dN_dx(jj));
               Klocal(TIp2, TJp1) += (b8*dN_dy(jj));
               Klocal(TIp2, TJp2) += 0.0;//(eps*b4*N(jj)); //

               // SUPG and PSPG stabilisation terms
               fact2 -= mu*d2N(jj);
               
               Dj(0,0) = FN(0,0) + fact2;
               Dj(0,1) = FN(0,1);
               Dj(0,2) = dN_dx(jj);
               Dj(1,0) = FN(1,0);
               Dj(1,1) = FN(1,1) + fact2;
               Dj(1,2) = dN_dy(jj);

               Klocal(TI, TJ)     += Da*Dj(0,0);
               Klocal(TI, TJp1)   += Da*Dj(0,1);
               Klocal(TI, TJp2)   += Da*Dj(0,2);

               Klocal(TIp1, TJ)   += Da*Dj(1,0);
               Klocal(TIp1, TJp1) += Da*Dj(1,1);
               Klocal(TIp1, TJp2) += Da*Dj(1,2);

              Klocal(TI,   TJ)   += ( (tau[0]*af*rho) * b1 * rStab(0) * N(jj) );  
              Klocal(TI,   TJp1) += ( (tau[0]*af*rho) * b2 * rStab(0) * N(jj) );
              Klocal(TIp1, TJ)   += ( (tau[0]*af*rho) * b1 * rStab(1) * N(jj) );
              Klocal(TIp1, TJp1) += ( (tau[0]*af*rho) * b2 * rStab(1) * N(jj) );

               Klocal(TIp2, TJ)   += (b1*Dj(0,0) + b2*Dj(1,0))*tau[1];
               Klocal(TIp2, TJp1) += (b1*Dj(0,1) + b2*Dj(1,1))*tau[1];
               Klocal(TIp2, TJp2) += (b1*Dj(0,2) + b2*Dj(1,2))*tau[1];
               
               // LSIC stabilisation
               
               //Klocal(TI,   TJ)   += (b1*dN_dx(jj))*tau[2];
               //Klocal(TI,   TJp1) += (b1*dN_dy(jj))*tau[2];

               //Klocal(TIp1, TJ)   += (b2*dN_dx(jj))*tau[2];
               //Klocal(TIp1, TJp1) += (b2*dN_dy(jj))*tau[2];

               if(axsy)
               {
                  // diffusion term
                  Klocal(TI, TJ)     += (mu*b8*N(jj)/rad/rad);
                  Klocal(TI, TJp2)   -= (b8*N(jj)/rad);

                  // continuity equation
                  Klocal(TIp2, TJ)   -= (b8*N(jj)/rad);
               }
             }

             Flocal(TI)   -= (b4*res2(0) + b1*stress(0,0) + b2*stress(0,1) );
             Flocal(TIp1) -= (b4*res2(1) + b1*stress(1,0) + b2*stress(1,1) );
             Flocal(TIp2) -= (b4*F.trace());

             // SUPG stabilisation terms
             Flocal(TI)   -= Da*rStab(0);
             Flocal(TIp1) -= Da*rStab(1);
             
             // PSPG stabilisation terms
             Flocal(TIp2) -= (tau[1]*(b1*rStab(0)+b2*rStab(1)));

             // LSIC stabilisation terms
             //Flocal(TI)   -= tau[2]*b1*F.trace();
             //Flocal(TIp1) -= tau[2]*b2*F.trace();

             if(axsy)
             {
                Flocal(TI)   += (-b4*(mu*vel(0)/rad/rad));
                Flocal(TI)   += (b4*pres/rad);
                Flocal(TIp2) += (b4*vel(0)/rad);
             }
          }
    }//gp1
    }//gp2
    
    return;
}
*/



template<>
void TreeNode<2>::applyBoundaryConditionsAtApoint(myDataIntegrateCutFEM& myData)
{
    int ii, jj, TI, nU;
   
    double  fact, af;
    VectorXd  NN(totnlbf), N;

    af = SolnData->td(2);

    GeomData->computeBasisFunctions2D(knotBegin, knotIncr, myData.param, NN);

    if(parent == NULL)
      N = NN;
    else
      N = SubDivMat*NN;
    
    //cout << " PENALTY " << PENALTY << endl;
    //Kovasznay analy;
    //myData.specVal[0] = analy.computeValue(2, myData.geom[0], myData.geom[1]);

    if(myData.dir < 3)
    {
      myData.specVal[0] -= computeValueCur(myData.dir, N);

      for(ii=0;ii<totnlbf2;ii++)
      {
        fact  = N[ii] * myData.dvol;

        TI = ndof*ii+myData.dir;

        myData.F1(TI) += (fact*myData.specVal[0]);
        fact *= af;

        for(jj=0;jj<totnlbf2;jj++)
          myData.K1(TI, ndof*jj+myData.dir) += fact * N[jj];
      }
    }
    if(myData.dir == 3)
    {
      nU  = forAssyVec.size();

      myData.specVal[0] -= computeValue2(0, N);

      for(ii=0;ii<totnlbf2;ii++)
      {
        fact = N[ii] * myData.dvol;

        TI = nU+ii;

        myData.F1(TI) += (fact*myData.specVal[0]);
        //fact *= af;

        for(jj=0;jj<totnlbf2;jj++)
        {
          myData.K1(TI, nU+jj) += fact * N[jj];
        }
      }
    }

    return;
}



                  /*
                  if(applied)
                  {
                       if(side == 6 & dir == 0)
                       {
                         if(uu<0.5)
                           val = tanh(50.0*uu);
                         else
                           val = -tanh(50.0*(uu-1.0));
                       }
                       else
                         val = DirichletData[side][dir];
                  }
                  else
                  {
                    if(dir == 0)
                      val = analy.computeXVelocity(uu, vv);
                    if(dir == 1)
                      val = analy.computeYVelocity(uu, vv);
                  }
                  if(side == 22)
                  {
                    if(dir == 1)
                      val = xx*(2.0-xx);
                    else
                      val = DirichletData[side][dir];
                  }

            if(side == 40)
            {
              if(dir == 0)
              {
                specVal = 10.0*(yy-0.25);
              }
            }
            if(side == 0)
            {
              if(dir == 0)
              {
                if(yy < y0 || yy > y1)
                  specVal = 0.0;
                else
                  specVal = specVal*(y1-yy)*(yy-y0);
              }
            }
            if(side == 33)
            {
              if(dir == 0)
              {
                //val = sin(PI*xx);
                if(uu<0.5)
                  specVal = tanh(50.0*uu);
                else
                  specVal = -tanh(50.0*(uu-1.0));
              }
            }
            if(side == 22)
            {
              if(dir == 1)
                specVal = xx*(2.0-xx);
            }
                  */


template<>
double TreeNode<1>::getJacBoundary(int side)
{
  double val;

  switch(side)
  {
      case 0:
            val = 0.5*knots[1][2] * GeomData->getJacobian(1);
      break;

      case 1:
            val = 0.5*knots[1][2] * GeomData->getJacobian(1);
      break;

      default :

            cout << " Invalid 'side' value in TreeNode<2>::applyBoundaryConditions " << endl;
      break;
  }

  return val;
}

template<>
double TreeNode<2>::getJacBoundary(int side)
{
  double val;

  switch(side)
  {
      case 0:
            val = 0.5*knots[1][2] * GeomData->getJacobian(1);
      break;

      case 1:
            val = 0.5*knots[1][2] * GeomData->getJacobian(1);
      break;

      case 2:
            val = 0.5*knots[0][2] * GeomData->getJacobian(0);
      break;

      case 3:
            val = 0.5*knots[0][2] * GeomData->getJacobian(0);
      break;

      default :

            cout << " Invalid 'side' value in TreeNode<2>::applyBoundaryConditions " << endl;
      break;
  }

  return val;
}


//
template<>
void TreeNode<2>::applyDirichletBCsGFEM(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
          //cout << "  ooooooooooooooo  " << endl;
  if( DirichletData.size() > 0 )
  {
    int ii, jj, aa, gp1, gp2, TI, TIp1, TIp2, index, TJ, TJp1, TJp2, dir, side, nU, nP;
    double theta, y0, y1, Ta, Tb, res, JacMultLoc, af, temp, NitscheFact;
    double  uu, vv, dvol, specVal, PENALTY, Jac, fact, rad, R, bb1, bb2;
    bool  isNitsche;

    myPoint  param, geom, normal, trac;
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
    //y0 = 0.0;    y1 = 0.5;

    //y0 = 0.1;    y1 = 0.6;
    //theta = 0.291456794;

    y0 = 0.0;    y1 = 1.61;
    //y0 = 0.0;    y1 = 0.41;
    //y0 = 0.0;    y1 = 2.0;
    //y0 = -5.0;    y1 = 5.0;

    R = 1.5;

    //Kovasznay  analy;
    //analy.SetPressure(1.310741966654558776639305506250821053981781005859375);
    //analy.SetPressure(0.0);

    //Stokes2DEx1  analy;
    //Stokes2DEx2  analy;
    //Stokes2DEx3  analy;


    //TwoDim_Ex1  analy;
    //PoissonEx2 analy;
    //PoissonEx1 analy;
    PoissonEx3 analy;
    //BiharmonicEx1 analy;

    //PearsonVortex analy;

      for(aa=0;aa<DirichletData.size();aa++)
      {
       // printVector(DirichletData[aa]);

        isNitsche = false;
        side        = (int) (DirichletData[aa][0] - 1);
        dir         = (int) (DirichletData[aa][1] - 1);
        specVal     = DirichletData[aa][2];
        PENALTY     = DirichletData[aa][3];
        isNitsche   = ( (int) DirichletData[aa][4] == 1 );
        NitscheFact = DirichletData[aa][5];

        //for symmetric Nitsche method -> NitscheFact = 1.0
        //for unsymmetric Nitsche method -> NitscheFact = -1.0

        //cout << " PENALTY " << PENALTY << endl;

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

            //R  = GeomData->computeCoord(0, 1.0);
            R = 0.5;
            geom[0] = GeomData->computeCoord(0, param[0]);
            geom[1] = GeomData->computeCoord(1, param[1]);
            //r = sqrt(xx*xx+yy*yy);
            //val = 1.0+log(2.0*r);
            //cout << xx << '\t' << yy << '\t' << DirichletData[aa][2] << endl;

            //cout << aa << '\t' << side << '\t' << dir << '\t' << val << endl;
            
            //res = specVal;

            if(axsy)
            {
              if(side > 0)
                dvol *= (2.0*PI*geom[0]);
            }

            specVal = DirichletData[aa][2];

            //
            if(side == 0 )
            {
              if(dir == 0)
              {
                //y0=0.0; y1=1.61;
                //specVal = DirichletData[aa][2]*(6.0/y1/y1)*(y1-geom[1])*(geom[1]-y0);
                //specVal = DirichletData[aa][2]*(y1-geom[1])*(geom[1]-y0);
                //cout << " specVal =  " << specVal << endl;
                
                y0=0.0; y1=0.5;
                //if(geom[1] <= 0.5)
                  specVal = DirichletData[aa][2]*(1.0-geom[1]*geom[1]/y1/y1);
                //else
                  //specVal = 0.0;
              }
            }
            //
            //
            //if(side == 3333)
            //{
              //if(dir == 0)
              //{
                //val = sin(PI*xx);
                //if(uu<0.5)
                  //specVal = tanh(20.0*uu);
                //else
                  //specVal = -tanh(20.0*(uu-1.0));
              //}
            //}
            //
            //
            //if(side == 1110 || side == 1111)
            //{
              //if(dir == 0)
              //{
                //specVal = DirichletData[aa][2]*(yy-0.25);
                //specVal = 1.5*yy*(2.0-yy);
                //if(yy <= y0 || yy >= y1)
                //if(xx > 1.5)
                  //specVal = 0.0;
                //else
                  //specVal = DirichletData[aa][2]*(y1-yy)*(yy-y0);
                  //specVal = (y1-xx)*(xx-y0);
                  //specVal = 3600.0*(1.0-yy*yy/R/R);
              //}
            //}
            //
            
            //if(side == 2222)
            //{
              //R = 0.5;
              //if(dir == 1)
              //{
                //if(xx >= R)
                  //specVal = 0.0;
                //else
                  //specVal = 1.0*xx*(1.5-xx);
                  //specVal = 3600.0*(1.0-xx*xx/R/R);
              //}
            //}
            //
            //

            //specVal = analy.computeValue(dir, xx, yy);
            //cout << side << '\t' << dir << '\t' << specVal << endl;

            //cout << " lllllllllll " << endl;
            specVal *= timeFunction[0].prop;
            //cout << " lllllllllll " << endl;
            //res *= mpapTime.cur;
            //specVal *= tanh(2.0*mpapTime.cur);
            //res *= sin(2.0*PI*mpapTime.cur+PI/2.0);
            //res *= sin(2.0*PI*10.0*mpapTime.cur+PI/2.0);
            //specVal *= 0.5*( 1.0-cos(628.3185*mpapTime.cur));

            //cout << aa << '\t' << side << '\t' << dir << '\t' << res << endl;
            //cout << specVal << '\t' << computeValue(dir, N) << endl;
            //cout << " PENALTY = " << PENALTY << endl;

              res = specVal - computeValue(dir, N);

              for(ii=0;ii<totnlbf2;ii++)
              {
                fact = N[ii] * dvol * PENALTY;

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
              res = specVal - computeValue(dir, N);

              //trac[0] = mu*(normal[0]*computeValue(0, dN_dx) + normal[1]*computeValue(0, dN_dy)) - normal[0]*computeValue2(0, N);
              trac[0] = mu*(normal[0]*computeValue(0, dN_dx) + normal[1]*computeValue(0, dN_dy)) - normal[0]*computeValue(2, N);
              for(ii=0;ii<totnlbf2;ii++)
              {
                TI   = ndof*ii;
                TIp1 = TI+1;
                TIp2 = TI+2;
                //TIp2 = nU+ii;

                Ta = mu*(normal[0]*dN_dx(ii)+normal[1]*dN_dy(ii))*dvol;
                bb1 = N[ii]*dvol;

                for(jj=0;jj<totnlbf2;jj++)
                {
                  TJ   = ndof*jj;
                  TJp1 = TJ+1;
                  TJp2 = TJ+2;
                  //TJp2 = nU+jj;

                  Tb = mu*(normal[0]*dN_dx(jj)+normal[1]*dN_dy(jj));
                  
                  fact = bb1*(-normal[0]*N(jj));

                  Klocal(TI, TJ)   -= (bb1*Tb);
                  Klocal(TI, TJp2) -= fact;

                  Klocal(TI, TJ)   -= (Ta*N(jj))*NitscheFact;
                  Klocal(TIp2, TJ) -= fact*NitscheFact;
                }

                Flocal(TI)   -= (-bb1*trac[0]);
                Flocal(TIp1) -= (0.0);

                Flocal(TI)   -= (Ta*res)*NitscheFact;
                Flocal(TIp2) -= (bb1*(-normal[0]*res))*NitscheFact;
              }
            } // if(dir == 0)
            if(dir == 1)
            {
              res = specVal - computeValue(1, N);

              //trac[1] = mu*(normal[0]*computeValue(1, dN_dx) + normal[1]*computeValue(1, dN_dy)) - normal[1]*computeValue2(0, N);
              trac[1] = mu*(normal[0]*computeValue(1, dN_dx) + normal[1]*computeValue(1, dN_dy)) - normal[1]*computeValue(2, N);
              for(ii=0;ii<totnlbf2;ii++)
              {
                TI   = ndof*ii;
                TIp1 = TI+1;
                TIp2 = TI+2;
                //TIp2 = nU+ii;

                Ta = mu*(normal[0]*dN_dx(ii)+normal[1]*dN_dy(ii))*dvol;
                bb1 = N[ii]*dvol;

                for(jj=0;jj<totnlbf2;jj++)
                {
                  TJ   = ndof*jj;
                  TJp1 = TJ+1;
                  TJp2 = TJ+2;
                  //TJp2 = nU+jj;

                  Tb = mu*(normal[0]*dN_dx(jj)+normal[1]*dN_dy(jj));
                  
                  fact = bb1*(-normal[1]*N(jj));

                  Klocal(TIp1, TJp1)  -= (bb1*Tb);
                  Klocal(TIp1, TJp2)  -= fact;

                  Klocal(TIp1, TJp1)  -= (Ta*N(jj))*NitscheFact;
                  Klocal(TIp2, TJp1)  -= fact*NitscheFact;
                }

                Flocal(TI)   -= (0.0);
                Flocal(TIp1) -= (-bb1*trac[1]);

                Flocal(TIp1) -= (Ta*res)*NitscheFact;
                Flocal(TIp2) -= (bb1*(-normal[1]*res))*NitscheFact;
              }
            } // if(dir == 1)
          } //if(isNitsche)

        }// for(gp1=0...
        }// for(gp2=0...
      } // for(aa=0;aa<DirichletData.size();aa++)
  } // if(DirichletData.size() > 0)

  return;
}
//




template<>
void TreeNode<2>::applyNeumannBCsGFEM(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
  if( NeumannData.size() > 0 )
  {
      int ii, jj, aa, gp1, gp2, TI, TIp1, TIp2, index, dir, side;
      double  y0, y1, res, JacMultLoc, af;
      double  dvol, specVal, PENALTY, Jac, rad, freq;

      VectorXd  N(totnlbf), dN_dx(totnlbf), dN_dy(totnlbf);
      VectorXd  NN(totnlbf), dNN_dx(totnlbf), dNN_dy(totnlbf);
      myPoint  param, geom;
      vector<double>  boundaryGPs1, boundaryGWs1, boundaryGPs2, boundaryGWs2;

      bool   axsy = ((int)elmDat[2] == 1);
      double  rho = elmDat[3];
      double  mu  = elmDat[4];

      af = SolnData->td(2);

      PENALTY  = 1.0;

      for(aa=0;aa<NeumannData.size();aa++)
      {
        side    = (int) (NeumannData[aa][0] - 1);
        dir     = (int) (NeumannData[aa][1] - 1);

        GeomData->setBoundaryGPs2D(side, boundaryGPs1, boundaryGWs1, boundaryGPs2, boundaryGWs2);

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
            }
            else
            {
              N = SubDivMat*NN;
            }

            //for(ii=0;ii<totnlbf;ii++)
            //printf(" \t %14.8f \n", NN[ii]);

            //printf("\n\n tracX and tracY ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n\n\n", xx, yy, val, val, Jac, JacMult, dvol);
            //if(pout) printf(" knotsAtGPs = %12.6f \t xx = %12.6f \t yy = %12.6f \n", knotsAtGPs, xx, yy);
            //printf(" uu = %12.6f \t vv = %12.6f \t dvol = %12.6f \t volume = %12.6f \n", uu, vv, dvol, volume);

            geom[0] = GeomData->computeCoord(0, param[0]);
            geom[1] = GeomData->computeCoord(1, param[1]);

            if(axsy)
              dvol *= (2.0*PI*geom[0]);

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
            //
            //
            if(side == 0 )
            {
              if(dir == 0)
              {
                if(geom[1] <= 0.5)
                  specVal = NeumannData[aa][2];
                else
                  specVal = 0.0;
              }
            }
            */

            res = dvol * specVal * timeFunction[0].prop;
            //res = dvol * specVal * tanh(20.0*mpapTime.cur);
            //cout << specVal << '\t' << axsy << '\t' << dir << endl;

            //freq = 2.0*PI*10.0;
            //res = dvol* specVal * (0.5*(1.0-cos(freq*mpapTime.cur)));
            //res = dvol * specVal * sin(freq*mpapTime.cur);
            //res *= tanh(50.0*timeFunction[0].prop);

            //cout << " yy " << yy << '\t' << res << endl;
            for(ii=0;ii<totnlbf2;ii++)
              Flocal(ndof*ii+dir) += (res*N(ii));

        }// for(gp1=0...
        }// for(gp2=0...
      } // for(aa=0;aa<NeumannData.size();aa++)
  } // if(NeumannData.size() > 0)

  return;
}


            /*
            if(dir == 0)
            {
              res = specVal - computeValue(dir, N);

              //trac[0] = mu*(normal[0]*computeValue(0, dN_dx) + normal[1]*computeValue(0, dN_dy)) - normal[0]*computeValue2(0, N);
              trac[0] = mu*(normal[0]*computeValue(0, dN_dx) + normal[1]*computeValue(0, dN_dy)) - normal[0]*computeValue(2, N);
              for(ii=0;ii<totnlbf2;ii++)
              {
                TI   = ndof*ii;
                TIp1 = TI+1;
                TIp2 = TI+2;
                //TIp2 = nU+ii;

                Ta = mu*(normal[0]*dN_dx(ii)+normal[1]*dN_dy(ii));

                for(jj=0;jj<totnlbf2;jj++)
                {
                  TJ   = ndof*jj;
                  TJp1 = TJ+1;
                  TJp2 = TJ+2;
                  //TJp2 = nU+jj;

                  Tb = mu*(normal[0]*dN_dx(jj)+normal[1]*dN_dy(jj));

                  Klocal(TI, TJ)   -= ((Ta*N(jj) + N(ii)*Tb)*dvol);
                  Klocal(TI, TJp2) -= (-normal[0]*N(ii)*N(jj)*dvol);
                  Klocal(TIp2, TJ) -= (-normal[0]*N(ii)*N(jj)*dvol);
                }

                Flocal(TI)   -= (dvol*(Ta*res-N(ii)*trac[0]));
                Flocal(TIp1) -= (0.0);
                Flocal(TIp2) -= (dvol*(-normal[0]*N(ii)*res));
              }
            }
            if(dir == 1)
            {
              res = specVal - computeValue(dir, N);

              //trac[1] = mu*(normal[0]*computeValue(1, dN_dx) + normal[1]*computeValue(1, dN_dy)) - normal[1]*computeValue2(0, N);
              trac[1] = mu*(normal[0]*computeValue(1, dN_dx) + normal[1]*computeValue(1, dN_dy)) - normal[1]*computeValue(2, N);
              for(ii=0;ii<totnlbf2;ii++)
              {
                TI   = ndof*ii;
                TIp1 = TI+1;
                TIp2 = TI+2;
                //TIp2 = nU+ii;

                Ta = mu*(normal[0]*dN_dx(ii)+normal[1]*dN_dy(ii));

                for(jj=0;jj<totnlbf2;jj++)
                {
                  TJ   = ndof*jj;
                  TJp1 = TJ+1;
                  TJp2 = TJ+2;
                  //TJp2 = nU+jj;

                  Tb = mu*(normal[0]*dN_dx(jj)+normal[1]*dN_dy(jj));

                  Klocal(TIp1, TJp1) -= ((Ta*N(jj) + N(ii)*Tb)*dvol);
                  Klocal(TIp1, TJp2) -= (-normal[1]*N(ii)*N(jj)*dvol);
                  Klocal(TIp2, TJp1) -= (-normal[1]*N(ii)*N(jj)*dvol);
                }

                Flocal(TI)   -= (0.0);
                Flocal(TIp1) -= (dvol*(Ta*res-N(ii)*trac[1]));
                Flocal(TIp2) -= (dvol*(-normal[1]*N(ii)*res));
              }
            }
            */
	    
	    
	    

